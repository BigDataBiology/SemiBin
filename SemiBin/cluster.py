import os
import math
import shutil
import tempfile

from .utils import write_bins, cal_num_bins
from .fasta import fasta_iter

# This is the default in the igraph package
NR_INFOMAP_TRIALS = 10

def run_infomap1(g, edge_weights, vertex_weights, trials):
    return g.community_infomap(edge_weights=edge_weights, vertex_weights=vertex_weights, trials=trials)

def run_infomap(g, edge_weights, vertex_weights, num_process):
    '''Run infomap, using multiple processors (if available)'''
    if num_process == 1:
        return g.community_infomap(edge_weights=edge_weights, vertex_weights=vertex_weights, trials=NR_INFOMAP_TRIALS)
    import multiprocessing
    with multiprocessing.Pool(num_process) as p:
        rs = [p.apply_async(run_infomap1, (g, edge_weights, vertex_weights, 1))
                for _ in range(NR_INFOMAP_TRIALS)]
        p.close()
        p.join()
    rs = [r.get() for r in rs]
    best = rs[0]
    for r in rs[1:]:
        if r.codelength < best.codelength:
            best = r
    return best



def cal_kl(m, v, use_ne='auto'):
    # A naive implementation creates a lot of copies of what can
    # become large matrices
    import numpy as np

    m = np.clip(m, 1e-6, None)
    v = np.clip(v, 1.0, None)

    m1 = m.reshape(1, len(m))
    m2 = m.reshape(len(m), 1)

    v1 = v.reshape(1, len(v))
    v2 = v.reshape(len(v), 1)


    if use_ne != 'no':
        try:
            import numexpr as ne

            res = ne.evaluate(
                    '(log(v1) - log(v2))/2 + ( (m1 - m2)**2 + v2 ) / ( 2 * v1 ) - half',
                    {
                        'v1': v1,
                        'v2': v2,
                        'm1': m1,
                        'm2': m2,
                        # numexpr rules are that mixing with floats causes
                        # conversion to float64
                        # Note that this does not happen for integers
                        'half': np.float32(0.5),
                    })
            np.clip(res, 1e-6, 1-1e-6, out=res)
            return res
        except ImportError:
            if use_ne != 'auto':
                raise
    v_div = np.log(v1) - np.log(v2)
    v_div /= 2.0

    m_dif = m1 - m2
    m_dif **=2
    m_dif += v2
    m_dif /= 2 * v1

    v_div += m_dif
    v_div -= 0.5
    np.clip(v_div, 1e-6, 1-1e-6, out=v_div)
    return v_div



def cluster(model, data, device, max_edges, max_node, is_combined,
            logger, n_sample, out, contig_dict, binned_length, num_process, minfasta, recluster, random_seed, orf_finder = 'prodigal'):
    """
    Cluster contigs into bins
    max_edges: max edges of one contig considered in binning
    max_node: max percentage of contigs considered in binning
    """
    from igraph import Graph
    from sklearn.neighbors import kneighbors_graph
    from sklearn.cluster import KMeans
    from scipy import sparse
    import torch
    import numpy as np
    train_data = data.values
    train_data_input = train_data[:, 0:136] if not is_combined else train_data

    depth = data.values[:, 136:len(data.values[0])].astype(np.float32)
    namelist = data.index.tolist()
    name2ix = {n:i for i, n in enumerate(namelist)}
    num_contigs = train_data_input.shape[0]
    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()
        embedding_matrix = kneighbors_graph(
            embedding,
            n_neighbors=min(max_edges, train_data.shape[0] - 1),
            mode='distance',
            p=2,
            n_jobs=num_process)
        kmer_matrix = kneighbors_graph(
            train_data_input,
            n_neighbors=min(max_edges, train_data.shape[0] - 1),
            mode='distance',
            p=2,
            n_jobs=num_process)

        # We want to intersect the matrices, so we make kmer_matrix into a
        # matrix of 1.s and multiply
        kmer_matrix.eliminate_zeros()
        kmer_matrix.data.fill(1.)
        embedding_matrix = embedding_matrix.multiply(kmer_matrix)
        embedding_matrix.eliminate_zeros()
        del kmer_matrix

    np.clip(embedding_matrix.data, None, 1., out=embedding_matrix.data)
    embedding_matrix.data = 1 - embedding_matrix.data

    threshold = 0
    max_axis1 = embedding_matrix.max(axis=1).toarray()
    while threshold < 1:
        threshold += 0.05
        n_above = np.sum(max_axis1 > threshold)
        if round(n_above / num_contigs, 2) < max_node:
            break
    threshold -= 0.05

    embedding_matrix.data[embedding_matrix.data<= threshold] = 0
    embedding_matrix.eliminate_zeros()
    if not is_combined:
        logger.info('Calculating depth matrix.')
        kl_matrix = None
        for k in range(n_sample):
            kl = cal_kl(depth[:,2*k], depth[:, 2*k + 1])
            if kl_matrix is None:
                kl *= -1
                kl_matrix = kl
            else:
                kl_matrix -= kl
        kl_matrix += n_sample
        kl_matrix /= n_sample
        embedding_matrix = embedding_matrix.multiply(kl_matrix)

    embedding_matrix.data[embedding_matrix.data <= 1e-6] = 0
    X, Y, V = sparse.find(embedding_matrix)
    # Find above diagonal positions:
    above_diag = Y > X
    X = X[above_diag]
    Y = Y[above_diag]

    edges = [(x,y) for x,y in zip(X, Y)]
    edge_weights = V[above_diag]

    logger.info('Edges:{}'.format(len(edges)))

    g = Graph()
    g.add_vertices(np.arange(num_contigs))
    g.add_edges(edges)
    length_weight = np.array([len(contig_dict[name]) for name in namelist])
    result = run_infomap(g,
                edge_weights=edge_weights,
                vertex_weights=length_weight,
                num_process=num_process)
    contig_labels = np.zeros(shape=num_contigs, dtype=int)

    for i, r in enumerate(result):
        contig_labels[r] = i

    output_bin_path = os.path.join(out, 'output_bins')
    if os.path.exists(output_bin_path):
        shutil.rmtree(output_bin_path)
    os.makedirs(output_bin_path, exist_ok=True)

    bin_files = write_bins(namelist, contig_labels, output_bin_path, contig_dict,minfasta=minfasta)
    if recluster:
        if not is_combined:
            mean_index = [2 * temp for temp in range(n_sample)]
            depth_mean = depth[:, mean_index]
            abun_scale = np.ceil(depth_mean.mean(axis = 0) / 100) * 100
            depth_mean = depth_mean / abun_scale
            scaling = np.mean(np.abs(embedding)) / np.mean(depth_mean)
            base = 10
            weight = 2 * base * math.ceil(scaling / base)
            embedding_new = np.concatenate(
                (embedding, depth_mean * weight), axis=1)
        else:
            embedding_new = embedding

        logger.info('Reclustering.')

        output_recluster_bin_path = os.path.join(out, 'output_recluster_bins')
        if os.path.exists(output_recluster_bin_path):
            shutil.rmtree(output_recluster_bin_path)

        os.makedirs(output_recluster_bin_path, exist_ok=True)

        with tempfile.TemporaryDirectory() as tdir:
            cfasta = os.path.join(tdir, 'concatenated.fna')
            with open(cfasta, 'wt') as concat_out:
                for ix,f in enumerate(bin_files):
                    for h,seq in fasta_iter(f):
                        h = f'bin{ix:06}.{h}'
                        concat_out.write(f'>{h}\n{seq}\n')
            seeds = cal_num_bins(
                cfasta,
                binned_length,
                num_process,
                multi_mode=True,
                orf_finder=orf_finder)

        for ix,bin_path in enumerate(bin_files):
            # if there are no hits, the output will be naturally empty
            seed = seeds.get(f'bin{ix:06}', [])
            num_bin = len(seed)

            if num_bin > 1:
                contig_list = [h for h,_ in fasta_iter(bin_path)]
                contig_index = [name2ix[c] for c in contig_list]
                re_bin_features = embedding_new[contig_index]

                seed_index = [name2ix[s] for s in seed]
                length_weight = np.array(
                    [len(contig_dict[name]) for name in contig_list])
                seeds_embedding = embedding_new[seed_index]
                kmeans = KMeans(
                    n_clusters=num_bin,
                    init=seeds_embedding,
                    n_init=1,
                    random_state=random_seed)
                kmeans.fit(re_bin_features, sample_weight=length_weight)
                labels = kmeans.labels_
                write_bins(contig_list, labels, os.path.join(out, 'output_recluster_bins'), contig_dict,
                           recluster=True, origin_label=int(bin_path.split('.')[-2]),minfasta = minfasta)
            else:
                shutil.copy(bin_path, os.path.join(out, 'output_recluster_bins'))

    logger.info('Binning finish.')

