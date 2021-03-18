import torch
from sklearn.neighbors import kneighbors_graph
from igraph import Graph
import numpy as np
from .utils import cal_kl, write_bins, cal_num_bins
import os
import math
from Bio import SeqIO
from sklearn.cluster import KMeans
import shutil


def cluster(model, data, device, max_edges, max_node, is_combined,
            logger, n_sample, contig_length_dict, out, contig_dict, binned_short,num_process,minfasta):
    train_data = data.values
    if not is_combined:
        train_data_input = train_data[:, 0:136]
    else:
        train_data_input = train_data

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
            n_neighbors=max_edges,
            mode='distance',
            p=2,
            n_jobs=-1).toarray()
        kmer_matrix = kneighbors_graph(
            train_data_input,
            n_neighbors=max_edges,
            mode='distance',
            p=2,
            n_jobs=-1).toarray()
        embedding_matrix[kmer_matrix == 0] = 0

    embedding_matrix[embedding_matrix >= 1] = 1
    embedding_matrix[embedding_matrix == 0] = 1
    embedding_matrix = 1 - embedding_matrix

    threshold = 0.95

    while (threshold >= 0):
        temp_threshold = embedding_matrix.copy()
        temp_threshold[temp_threshold <= threshold] = 0
        num = len(list(set(np.where(temp_threshold > 0)[0])))
        if round(num / len(embedding_matrix), 2) >= max_node:
            break
        else:
            threshold -= 0.05

    embedding_matrix[embedding_matrix <= threshold] = 0
    if not is_combined:
        logger.info('Calculating depth matrix.')
        depth_matrix = np.zeros(shape=embedding_matrix.shape)
        for i in range(len(embedding_matrix)):
            for j in range(i + 1, len(embedding_matrix)):
                if embedding_matrix[i][j] > 0:
                    temp_depth = 0
                    for k in range(n_sample):
                        temp_depth += 1 - \
                            cal_kl(depth[i][2 * k], depth[j][2 * k],
                                   depth[i][2 * k + 1], depth[j][2 * k + 1])
                    depth_matrix[i][j] = temp_depth / n_sample

        matrix = embedding_matrix * depth_matrix
    else:
        matrix = embedding_matrix

    edges = []
    edges_weight = []

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if matrix[i][j] > 1e-6:
                edges.append((i, j))
                edges_weight.append(matrix[i][j])

    logger.info('Edges:{}'.format(len(edges)))

    g = Graph()
    g.add_vertices(list(range(len(matrix))))
    g.add_edges(edges)
    length_weight = np.array([contig_length_dict[name] for name in namelist])
    result = g.community_infomap(
        edge_weights=edges_weight,
        vertex_weights=length_weight)
    contig_labels = np.zeros(shape=(len(matrix)), dtype=np.int)

    for i, r in enumerate(result):
        for infomap_index in r:
            contig_labels[infomap_index] = i

    output_bin_path = os.path.join(out, 'output_bins')
    os.makedirs(output_bin_path, exist_ok=True)

    write_bins(namelist, contig_labels, output_bin_path, contig_dict,minfasta=minfasta)
    if not is_combined:
        mean_index = [2 * temp for temp in range(n_sample)]
        depth_mean = depth[:, mean_index] / 100
        scaling = np.mean(np.abs(embedding)) / np.mean(depth_mean)
        base = 10
        weight = 2 * base * math.ceil(scaling / base)
        embedding_new = np.concatenate(
            (embedding, depth_mean * weight), axis=1)
    else:
        embedding_new = embedding

    bin_files = os.listdir(output_bin_path)
    logger.info('Reclustering.')

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
                    binned_short,num_process)
            except BaseException:
                pass
            contig_index = [mapObj[temp] for temp in contig_list]
            re_bin_features = embedding_new[contig_index]
            if not os.path.exists(os.path.join(out, 'output_recluster_bins')):
                os.mkdir(os.path.join(out, 'output_recluster_bins'))

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

