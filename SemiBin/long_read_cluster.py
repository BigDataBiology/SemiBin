import os
import torch
import numpy as np
from .utils import cal_num_bins, get_marker, write_bins
from sklearn.cluster import DBSCAN
from sklearn.neighbors import kneighbors_graph
from collections import defaultdict


def get_max(results_dict, con, contig_to_marker, namelist, contig_dict, minfasta):
    max_F1 = 0
    max_weight = 1e9
    max_bin = None
    for eps_value, res_labels in results_dict.items():
        res = defaultdict(list)
        for label, name in zip(res_labels, namelist):
            if label != -1:
                res[label].append(name)
        for temp in res:
            bin_contig = res[temp]
            if sum(len(contig_dict[contig]) for contig in
                   bin_contig) < minfasta:
                continue
            marker_list = []
            for contig in bin_contig:
                marker_list.extend(contig_to_marker[contig])
            if len(marker_list) == 0:
                continue
            recall = len(set(marker_list)) / 107
            contamination = (len(marker_list) - len(set(marker_list))) / len(
                marker_list)
            if contamination <= con:
                F1 = 2 * recall * (1 - contamination) / (
                            recall + (1 - contamination))
                if F1 > max_F1:
                    max_F1 = F1
                    max_weight = sum(
                        len(contig_dict[contig]) for contig in bin_contig)
                    max_bin = bin_contig
                if F1 == max_F1:
                    if sum(len(contig_dict[contig]) for contig in
                           bin_contig) <= max_weight:
                        max_weight = sum(
                            len(contig_dict[contig]) for contig in bin_contig)
                        max_bin = bin_contig
    return max_F1, max_weight, max_bin

def cluster_long_read(model, data, device, is_combined,
            logger, n_sample, out, contig_dict, binned_length, num_process, minfasta, random_seed, orf_finder = 'prodigal'):
    contig_list = data.index.tolist()
    if not is_combined:
        train_data_input = data.values[:, 0:136]
    else:
        train_data_input = data.values
        if train_data_input.shape[1] - 136 > 20:
            train_data_kmer = train_data_input[:, 0:136]
            train_data_depth = train_data_input[:, 136:len(data.values[0])]
            from sklearn.preprocessing import normalize
            train_data_depth = normalize(train_data_depth, axis=1, norm='l1')
            train_data_input = np.concatenate((train_data_kmer, train_data_depth), axis=1)

    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()

    length_weight = np.array(
        [len(contig_dict[name]) for name in contig_list])

    if not is_combined:
        depth = data.values[:, 136:len(data.values[0])].astype(np.float32)
        mean_index = [2 * temp for temp in range(n_sample)]
        depth = depth[:, mean_index]
        embedding_new = np.concatenate((embedding, np.log(depth)), axis=1)
    else:
        embedding_new = embedding

    cfasta = os.path.join(out, 'concatenated.fna')
    with open(cfasta, 'wt') as concat_out:
        for h in contig_list:
            concat_out.write(f'>{h}\n{contig_dict[h]}\n')
    seeds = cal_num_bins(
        cfasta,
        binned_length,
        64,
        output=out,
        orf_finder=orf_finder)

    contig2marker = get_marker(f'{out}/markers.hmmout', orf_finder=orf_finder,
                               min_contig_len=binned_length, fasta_path=cfasta, contig_to_marker=True)

    output_bin_path = os.path.join(out, 'output_bins')

    logger.debug('Running DBSCAN.')
    dist_matrix = kneighbors_graph(
        embedding_new,
        n_neighbors=min(200, embedding_new.shape[0] - 1),
        mode='distance',
        p=2,
        n_jobs=num_process)

    DBSCAN_results_dict = {}
    for eps_value in [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55]:
        dbscan = DBSCAN(eps=eps_value, min_samples=5, n_jobs=num_process, metric='precomputed')
        dbscan.fit(dist_matrix, sample_weight=length_weight)
        labels = dbscan.labels_
        DBSCAN_results_dict[eps_value] = labels.tolist()

    cluster_label = 0
    logger.debug('Integrating results.')

    while sum(len(contig_dict[contig]) for contig in contig_list) >= minfasta:
        if len(contig_list) == 1:
            write_bins(contig_list,
                       [cluster_label] * len(contig_list),
                       output_bin_path, contig_dict,
                       recluster=False,
                       minfasta=minfasta)
            cluster_label += 1
            break

        for con in [0.1, 0.2, 0.3, 0.4, 0.5, 1]:
            max_F1, max_weight, max_bin = get_max(DBSCAN_results_dict, con, contig2marker, contig_list, contig_dict, minfasta)
            if max_F1 != 0:
                break
            else:
                if con != 1:
                    continue
                else:
                    break

        if max_bin == [] or max_bin == None:
            break
        else:
            write_bins(max_bin, [cluster_label] * len(max_bin),
                       output_bin_path, contig_dict,
                       recluster=False,
                       minfasta=minfasta)
            cluster_label += 1

            for temp in max_bin:
                temp_index = contig_list.index(temp)
                contig_list.pop(temp_index)
                for eps_value in DBSCAN_results_dict:
                    DBSCAN_results_dict[eps_value].pop(temp_index)

    logger.info('Finished binning.')