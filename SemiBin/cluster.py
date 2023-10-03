import os
from os import path
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



def run_embed_infomap(logger, model, data, * ,
            device, max_edges, max_node, is_combined, n_sample, contig_dict,
            num_process : int, random_seed):
    """
    Cluster contigs into bins
    max_edges: max edges of one contig considered in binning
    max_node: max percentage of contigs considered in binning
    """
    from igraph import Graph
    from sklearn.neighbors import kneighbors_graph
    from scipy import sparse
    import torch
    import numpy as np

    train_data_input = data.values[:, 0:136] if not is_combined else data.values

    if is_combined:
        if train_data_input.shape[1] - 136 > 20:
            train_data_kmer = train_data_input[:, 0:136]
            train_data_depth = train_data_input[:, 136:len(data.values[0])]
            from sklearn.preprocessing import normalize
            train_data_depth = normalize(train_data_depth, axis=1, norm='l1')
            train_data_input = np.concatenate((train_data_kmer, train_data_depth), axis=1)

    depth = data.values[:, 136:len(data.values[0])].astype(np.float32)
    num_contigs = train_data_input.shape[0]
    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()
        embedding_matrix = kneighbors_graph(
            embedding,
            n_neighbors=min(max_edges, data.shape[0] - 1),
            mode='distance',
            p=2,
            n_jobs=num_process)
        kmer_matrix = kneighbors_graph(
            train_data_input,
            n_neighbors=min(max_edges, data.shape[0] - 1),
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
        logger.debug('Calculating depth matrix.')
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

    logger.debug(f'Number of edges in clustering graph: {len(edges)}')

    g = Graph()
    g.add_vertices(np.arange(num_contigs))
    g.add_edges(edges)
    length_weight = np.array([len(contig_dict[name]) for name in data.index])
    logger.debug(f'Running infomap with {num_process} processes...')
    result = run_infomap(g,
                edge_weights=edge_weights,
                vertex_weights=length_weight,
                num_process=num_process)
    contig_labels = np.zeros(shape=num_contigs, dtype=int)

    for i, r in enumerate(result):
        contig_labels[r] = i
    return embedding, contig_labels


def recluster_bins(logger, data, *, n_sample, embedding, is_combined, contig_labels, minfasta, contig_dict, binned_length, orf_finder, num_process, random_seed):
    from sklearn.cluster import KMeans
    import numpy as np
    from collections import defaultdict
    if not is_combined:
        depth = data.values[:, 136:len(data.values[0])].astype(np.float32)
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


    total_size = defaultdict(int)
    for i, c in enumerate(contig_labels):
        total_size[c] += len(contig_dict[data.index[i]])
    with tempfile.TemporaryDirectory() as tdir:
        cfasta = os.path.join(tdir, 'concatenated.fna')
        with open(cfasta, 'wt') as concat_out:
            for ix,h in enumerate(data.index):
                bin_ix = contig_labels[ix]
                if total_size[bin_ix] < minfasta:
                    continue
                concat_out.write(f'>bin{bin_ix:06}.{h}\n')
                concat_out.write(contig_dict[data.index[ix]] + '\n')

        seeds = cal_num_bins(
            cfasta,
            binned_length,
            num_process,
            multi_mode=True,
            orf_finder=orf_finder)
            # we cannot bypass the orf_finder here, because of the renaming of the contigs
        if seeds == []:
            logger.warning('No bins found in the concatenated fasta file.')
            return contig_labels

    name2ix = {name:ix for ix,name in enumerate(data.index)}
    contig_labels_reclustered = np.empty_like(contig_labels)
    contig_labels_reclustered.fill(-1)
    next_label = 0
    for bin_ix in range(contig_labels.max() + 1):
        seed = seeds.get(f'bin{bin_ix:06}', [])
        num_bin = len(seed)

        if num_bin > 1 and total_size[bin_ix] >= minfasta:
            contig_indices = [i for i,ell in enumerate(contig_labels) if ell == bin_ix]
            re_bin_features = embedding_new[contig_indices]

            seed_index = [name2ix[s] for s in seed]
            length_weight = np.array(
                [len(contig_dict[name])
                    for name,ell in zip(data.index, contig_labels)
                        if ell == bin_ix])
            seeds_embedding = embedding_new[seed_index]
            kmeans = KMeans(
                n_clusters=num_bin,
                init=seeds_embedding,
                n_init=1,
                random_state=random_seed)
            kmeans.fit(re_bin_features, sample_weight=length_weight)
            for i, label in enumerate(kmeans.labels_):
                contig_labels_reclustered[contig_indices[i]] = next_label + label
            next_label += num_bin
        else:
            contig_labels_reclustered[contig_labels == bin_ix] = next_label
            next_label += 1
    assert contig_labels_reclustered.min() >= 0
    return contig_labels_reclustered


def cluster(logger, model, data, device, is_combined,
            n_sample, out, contig_dict, *, args,
            binned_length, minfasta):
    """
    Cluster contigs into bins
    max_edges: max edges of one contig considered in binning
    max_node: max percentage of contigs considered in binning
    """
    import pandas as pd
    import numpy as np
    embedding, contig_labels = run_embed_infomap(logger, model, data,
            device=device, max_edges=args.max_edges, max_node=args.max_node,
            is_combined=is_combined, n_sample=n_sample,
            contig_dict=contig_dict, num_process=args.num_process,
            random_seed=args.random_seed)

    if args.write_pre_reclustering_bins or not args.recluster:
        output_bin_path = os.path.join(out,
                    'output_prerecluster_bins' if args.recluster else 'output_bins')
        if os.path.exists(output_bin_path):
            logger.warning(f'Previous output directory `{output_bin_path}` found. Over-writing it.')
            shutil.rmtree(output_bin_path)
        os.makedirs(output_bin_path, exist_ok=True)


        bin_files = write_bins(data.index.tolist(),
                            contig_labels,
                            output_bin_path,
                            contig_dict,
                            minfasta=minfasta,
                            output_tag=args.output_tag,
                            output_compression=args.output_compression)
        if not len(bin_files):
            logger.warning('No bins were created. Please check your input data.')
            return
        if not args.recluster:
            logger.info(f'Number of bins: {len(bin_files)}')
            bin_files.to_csv(os.path.join(out, 'bins_info.tsv'), index=False, sep='\t')
            pd.DataFrame({'contig': data.index, 'bin': contig_labels}).to_csv(
                os.path.join(out, 'contig_bins.tsv'), index=False, sep='\t')
        n_pre_bins = len(bin_files)
    else:
        from collections import defaultdict
        total_size = defaultdict(int)
        for i, c in enumerate(contig_labels):
            total_size[c] += len(contig_dict[data.index[i]])
        n_pre_bins = sum((total_size[bin_ix] >= minfasta for bin_ix in range(contig_labels.max() + 1)))


    if args.recluster:
        logger.info(f'Number of bins prior to reclustering: {n_pre_bins}')
        logger.debug('Reclustering...')

        contig_labels_reclustered = recluster_bins(logger,
                                                data,
                                                n_sample=n_sample,
                                                embedding=embedding,
                                                contig_labels=contig_labels,
                                                contig_dict=contig_dict,
                                                minfasta=minfasta,
                                                binned_length=binned_length,
                                                num_process=args.num_process,
                                                orf_finder=args.orf_finder,
                                                random_seed=args.random_seed,
                                                is_combined=is_combined)
        output_recluster_bin_path = path.join(out,
                        ('output_recluster_bins'
                            if args.write_pre_reclustering_bins
                            else 'output_bins'))
        if os.path.exists(output_recluster_bin_path):
            logger.warning(f'Previous output directory `{output_recluster_bin_path}` found. Over-writing it.')
            shutil.rmtree(output_recluster_bin_path)
        os.makedirs(output_recluster_bin_path, exist_ok=True)
        outputs = write_bins(data.index.tolist(),
                            contig_labels_reclustered,
                            output_recluster_bin_path,
                            contig_dict,
                            minfasta=minfasta,
                            output_tag=args.output_tag,
                            output_compression=args.output_compression)
        logger.info(f'Number of bins after reclustering: {len(outputs)}')
        outputs.to_csv(os.path.join(out, 'recluster_bins_info.tsv'), index=False, sep='\t')
        pd.DataFrame({'contig': data.index, 'bin': contig_labels_reclustered}).to_csv(
            os.path.join(out, 'contig_bins.tsv'), index=False, sep='\t')
    logger.info('Binning finished')


