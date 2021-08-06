from sklearn.neighbors import kneighbors_graph
from igraph import Graph
from .utils import cal_kl, write_bins, cal_num_bins
import os
import math
from Bio import SeqIO
from sklearn.cluster import KMeans
import shutil


def cluster(model, data, device, max_edges, max_node, is_combined,
            logger, n_sample, contig_length_dict, out, contig_dict, binned_length,num_process,minfasta,recluster,random_seed):
    """
    Cluster contigs into bins
    max_edges: max edges of one contig considered in binning
    max_node: max percentage of contigs considered in binning
    """
    import torch
    import numpy as np
    train_data = data.values
    train_data_input = train_data[:, 0:136] if not is_combined else train_data

    depth = data.values[:, 136:len(data.values[0])]
    namelist = data.index.tolist()
    row_index = data._stat_axis.values.tolist()
    mapObj = {n:i for i, n in enumerate(namelist)}
    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()
        embedding_matrix = kneighbors_graph(
            embedding,
            n_neighbors=min(max_edges, train_data.shape[0] - 1),
            mode='distance',
            p=2,
            n_jobs=num_process).toarray()
        kmer_matrix = kneighbors_graph(
            train_data_input,
            n_neighbors=min(max_edges, train_data.shape[0] - 1),
            mode='distance',
            p=2,
            n_jobs=num_process).toarray()
        embedding_matrix[kmer_matrix == 0] = 0
    kmer_matrix = None
    del kmer_matrix

    embedding_matrix[embedding_matrix >= 1] = 1
    embedding_matrix[embedding_matrix == 0] = 1
    embedding_matrix = 1 - embedding_matrix

    threshold = 0

    while (threshold < 1):
        threshold += 0.05
        num = len(list(set(np.where(embedding_matrix > threshold)[0])))
        if round(num / len(embedding_matrix), 2) < max_node:
            break
    threshold -= 0.05

    embedding_matrix[embedding_matrix <= threshold] = 0
    if not is_combined:
        logger.info('Calculating depth matrix.')
        kl_matrix = np.zeros(shape=(len(embedding_matrix), len(embedding_matrix)))
        for k in range(n_sample):
            kl = cal_kl(depth[:,2*k], depth[:, 2*k + 1])
            kl_matrix +=  1 - kl
        kl_matrix = kl_matrix / n_sample
        embedding_matrix = embedding_matrix * kl_matrix

    edges = []
    edges_weight = []

    for i in range(len(embedding_matrix)):
        for j in range(i + 1, len(embedding_matrix)):
            if embedding_matrix[i][j] > 1e-6:
                edges.append((i, j))
                edges_weight.append(embedding_matrix[i][j])

    logger.info('Edges:{}'.format(len(edges)))

    g = Graph()
    g.add_vertices(list(range(len(embedding_matrix))))
    g.add_edges(edges)
    length_weight = np.array([contig_length_dict[name] for name in namelist])
    result = g.community_infomap(edge_weights=edges_weight, vertex_weights=length_weight)
    contig_labels = np.zeros(shape=(len(embedding_matrix)), dtype=np.int)

    for i, r in enumerate(result):
        for infomap_index in r:
            contig_labels[infomap_index] = i

    output_bin_path = os.path.join(out, 'output_bins')
    if os.path.exists(output_bin_path):
        shutil.rmtree(output_bin_path)
    os.makedirs(output_bin_path, exist_ok=True)

    write_bins(namelist, contig_labels, output_bin_path, contig_dict,minfasta=minfasta)
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

        bin_files = os.listdir(output_bin_path)
        logger.info('Reclustering.')

        output_recluster_bin_path = os.path.join(out, 'output_recluster_bins')
        if os.path.exists(output_recluster_bin_path):
            shutil.rmtree(output_recluster_bin_path)

        os.makedirs(output_recluster_bin_path, exist_ok=True)

        for bin in bin_files:
            if os.path.exists(os.path.join(output_bin_path, bin)):
                contig_list = []
                for seq_record in SeqIO.parse(
                        os.path.join(output_bin_path, bin), "fasta"):
                    contig_list.append(seq_record.id)
                contig_output = os.path.join(output_bin_path, bin) + '.frag'
                hmm_output = os.path.join(output_bin_path, bin) + '.hmmout'
                seed_output = os.path.join(output_bin_path, bin) + '.seed'
                try:
                    cal_num_bins(
                        os.path.join(output_bin_path,bin),
                        contig_output,
                        hmm_output,
                        seed_output,
                        binned_length,
                        num_process)
                except BaseException:
                    pass
                contig_index = [mapObj[temp] for temp in contig_list]
                re_bin_features = embedding_new[contig_index]

                if os.path.exists(seed_output):
                    seed = open(seed_output).read().split('\n')
                    seed = [contig for contig in seed if contig != '']
                    init_seed = seed
                    num_bin = len(seed)
                    seed_index = []
                    for temp in init_seed:
                        seed_index.append(row_index.index(temp))
                    length_weight = np.array(
                        [contig_length_dict[name] for name in contig_list])
                    seeds_embedding = embedding_new[seed_index]
                    if random_seed is not None:
                        kmeans = KMeans(
                            n_clusters=num_bin,
                            init=seeds_embedding,
                            n_init=1,
                            random_state=random_seed)
                    else:
                        kmeans = KMeans(
                            n_clusters=num_bin,
                            init=seeds_embedding,
                            n_init=1)
                    kmeans.fit(re_bin_features, sample_weight=length_weight)
                    labels = kmeans.labels_
                    write_bins(contig_list, labels, os.path.join(out, 'output_recluster_bins'), contig_dict,
                               recluster=True, origin_label=int(bin.split('.')[-2]),minfasta = minfasta)
                else:
                    shutil.copy(os.path.join(
                        output_bin_path, bin), os.path.join(out, 'output_recluster_bins'))

    logger.info('Binning finish.')

