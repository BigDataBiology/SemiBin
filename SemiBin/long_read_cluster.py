import os
import torch
import numpy as np
import polars as pl
from .utils import cal_num_bins, get_marker, write_bins, normalize_kmer_motif_features
from sklearn.cluster import DBSCAN
from sklearn.neighbors import kneighbors_graph
from collections import defaultdict
from scipy.sparse import save_npz

def get_best_bin(results_dict, contig_to_marker, namelist, contig_dict, minfasta):

    # There is room for improving the loop below to avoid repeated computation
    # but it runs very fast in any case
    for max_contamination in [0.1, 0.2, 0.3, 0.4, 0.5, 1]:
        max_F1 = 0
        weight_of_max = 1e9
        max_bin = None

        for res_labels in results_dict.values():
            res = defaultdict(list)
            for label, name in zip(res_labels, namelist):
                if label != -1:
                    res[label].append(name)
            for bin_contig in res.values():
                cur_weight = sum(len(contig_dict[contig]) for contig in bin_contig)
                if cur_weight < minfasta:
                    continue
                marker_list = []
                for contig in bin_contig:
                    marker_list.extend(contig_to_marker[contig])
                if len(marker_list) == 0:
                    continue
                recall = len(set(marker_list)) / 107
                contamination = (len(marker_list) - len(set(marker_list))) / len(
                    marker_list)
                if contamination <= max_contamination:
                    F1 = 2 * recall * (1 - contamination) / (
                                recall + (1 - contamination))
                    if F1 > max_F1:
                        max_F1 = F1
                        weight_of_max = cur_weight
                        max_bin = bin_contig
                    elif F1 == max_F1 and cur_weight <= weight_of_max:
                        weight_of_max = cur_weight
                        max_bin = bin_contig
        if max_F1 > 0: # if there is a bin with F1 > 0
            return max_bin

# Remove clusters smaller than 500_000 bp
def remove_small_clusters(extracted, contig_dict, unbinned_df, results_df, min_fasta, logger):
    logger.info("Removing clusters smaller than {} bp".format(min_fasta))
    n_clusters = len(extracted)
    logger.info("Number of clusters before removal: {}".format(n_clusters))
    for cluster in extracted:
        if sum(len(contig_dict[contig]) for contig in cluster) < min_fasta:
            extracted.remove(cluster)
            eps_values_for_cluster = results_df.filter(pl.col("Contig").is_in(cluster))
            unbinned_df = pl.concat([unbinned_df, eps_values_for_cluster])
            assert unbinned_df.get_column("Contig").n_unique() == unbinned_df.shape[0], "Not all contigs are unique"    
    logger.info("removed {} clusters".format(n_clusters - len(extracted)))
    return extracted, unbinned_df



def cluster_long_read(logger, model, data, device, is_combined,
            n_sample, out, contig_dict, *, binned_length, args,
            minfasta, features_data):
    import pandas as pd
    from .utils import norm_abundance
    contig_list = data.index.tolist()
    if not is_combined:
        train_data_input = data.values[:, features_data["kmer"] + features_data["motif"]]
        train_data_input, _ = normalize_kmer_motif_features(train_data_input, train_data_input)
    else:
        train_data_input = data.values
        if norm_abundance(train_data_input, features_data):
            train_data_kmer = train_data_input[:, features_data["kmer"] + features_data["motif"]]
            train_data_kmer, _ = normalize_kmer_motif_features(train_data_kmer, train_data_kmer)
            
            train_data_depth = train_data_input[:, features_data["depth"]]
            from sklearn.preprocessing import normalize
            train_data_depth = normalize(train_data_depth, axis=1, norm='l1')
            train_data_input = np.concatenate((train_data_kmer, train_data_depth), axis=1)

    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()
        

    length_weight = np.array(
        [len(contig_dict[name]) for name in contig_list])
    length_weight = np.log10(length_weight)
        
    if not is_combined:
        depth = data.values[:, features_data["depth"]].astype(np.float32)
        mean_index = [2 * temp for temp in range(n_sample)]
        depth = depth[:, mean_index]
        embedding_new = np.concatenate((embedding, np.log(depth)), axis=1)
    else:
        embedding_new = embedding

    # Create a DataFrame to save the embeddings
    embedding_df = pd.DataFrame(embedding_new, index=contig_list)

    # Save the embeddings to a CSV file
    embedding_df.to_csv(os.path.join(out, 'embedding_space.csv'))

    import tempfile
    with tempfile.TemporaryDirectory() as tdir:
        cfasta = os.path.join(tdir, 'concatenated.fna')
        with open(cfasta, 'wt') as concat_out:
            for h in contig_list:
                concat_out.write(f'>{h}\n{contig_dict[h]}\n')
            seeds = cal_num_bins(
                cfasta,
                binned_length,
                args.num_process,
                output=out,
                orf_finder=args.orf_finder,
                prodigal_output_faa=args.prodigal_output_faa)

            contig2marker = get_marker(os.path.join(out, 'markers.hmmout'), orf_finder=args.orf_finder,
                                       min_contig_len=binned_length, fasta_path=cfasta, contig_to_marker=True)

    output_bin_path = os.path.join(out, 'output_bins')

    logger.debug('Running DBSCAN.')
    dist_matrix = kneighbors_graph(
        embedding_new,
        n_neighbors=min(200, embedding_new.shape[0] - 1),
        mode='distance',
        p=2,
        n_jobs=args.num_process)

    # Save distance matrix
    save_npz(os.path.join(out, "dist_matrix.npz"), dist_matrix)
    
    DBSCAN_results_dict = {}
    eps_values = [0.01] + [round(x, 2) for x in np.arange(0.05, 1, 0.05)]
    for eps_value in eps_values:
        dbscan = DBSCAN(eps=eps_value, min_samples=5, n_jobs=args.num_process, metric='precomputed')
        dbscan.fit(dist_matrix, sample_weight=length_weight)
        labels = dbscan.labels_
        DBSCAN_results_dict[eps_value] = labels.tolist()

    
    # Prepare DataFrame for saving results
    results_df = pd.DataFrame({'Contig': contig_list})

    # Add cluster labels for each eps value to the DataFrame
    for eps_value in eps_values:
        results_df[f'Cluster_Label_eps_{eps_value}'] = DBSCAN_results_dict[eps_value]

    # Save results to CSV
    results_df.to_csv(os.path.join(out,'dbscan_results_multiple_eps.csv'), index=False)
    results_df = pl.DataFrame(results_df)
    
    logger.debug('Integrating results.')

    extracted = []
    initial_cluster_dict = {k: v for k, v in DBSCAN_results_dict.items() if 0.01 <= k <= 0.55}
    # while the sum of contigs is higher than the smallest bin allowed continue to find bins
    while sum(len(contig_dict[contig]) for contig in contig_list) >= minfasta:
        # If there is only one contig left, add it to the extracted list
        if len(contig_list) == 1:
            extracted.append(contig_list)
            break
        
        
        max_bin = get_best_bin(initial_cluster_dict,
                                contig2marker,
                                contig_list,
                                contig_dict,
                                minfasta)
        if not max_bin:
            break

        extracted.append(max_bin)
        for temp in max_bin:
            temp_index = contig_list.index(temp)
            contig_list.pop(temp_index)
            for eps_value in DBSCAN_results_dict:
                DBSCAN_results_dict[eps_value].pop(temp_index)
        
        
    contig2ix = {}
    for i, cs in enumerate(extracted):
        for c in cs:
            contig2ix[c] = i
    namelist = data.index.tolist()
    contig_labels = [contig2ix.get(c, -1) for c in namelist]
    written = write_bins(namelist, contig_labels,
                    output_bin_path, contig_dict,
                    minfasta=minfasta,
                    output_tag=args.output_tag,
                    output_compression=args.output_compression)
    logger.info(f'Number of bins: {len(written)}')
    written.to_csv(os.path.join(out, 'bins_info.tsv'), index=False,
                   sep='\t')
    pd.DataFrame({'contig': namelist, 'bin': contig_labels}).to_csv(
            os.path.join(out, 'contig_bins.tsv'), index=False, sep='\t')
    logger.info('Finished binning.')
