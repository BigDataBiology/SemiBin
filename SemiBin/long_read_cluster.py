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



def expand_cluster_based_on_methylation(
    extracted,
    results_df,
    unbinned_df,
    eps_values,
    data,
    features_data,
    logger,
    args,
    
):
    from .methylation_pattern import generate_methylation_pattern
    # expand clusters
    motifs = data.columns.tolist()[min(features_data["motif"]):(max(features_data["motif"]) + 1)]
    
    
    methylation_comparison = generate_methylation_pattern(
        logger = logger,
        motifs_scored = args.motifs_scored,
        motif_features = motifs,
        ambiguous_interval = [0.1, 0.6],
        min_motif_observations = 5,
        threads = 20,
        min_comparisons = 5
    )
    
    refined_extracted = []
    t_contigs_added = 0
    for ibin in extracted:
        logger.debug("Processing ibin: {}".format(ibin))
        ibin_eps = results_df\
            .filter(pl.col("Contig").is_in(ibin))\
            .melt(id_vars = "Contig")\
            .filter(pl.col("value") != -1)\
            .group_by("variable")\
            .agg(
                pl.col("value").n_unique().alias("unique_values"),
                pl.col("value").len().alias("count")    
            )\
            .filter((pl.col("unique_values") == 1) & (pl.col("count") == len(ibin)))\
            .with_columns(
                pl.col("variable").str.extract(r"_(\d+\.\d+)").cast(pl.Float32).alias("eps")
            )\
            .sort("eps")\
            .get_column("variable")[0]
        # strip string from "Cluster_label_eps_" to get the eps value and convert to float
        ibin_eps = ibin_eps[18:]
        
        # Check methylation pattern of ibin
        ibin_pattern = methylation_comparison\
            .filter(pl.col("contig").is_in(ibin))\
            .filter(pl.col("contig_compare").is_in(ibin))
            
        # Check the ibin for mismatches
        if sum(ibin_pattern["sum_mismatches"] > 0):
            logger.debug("Mismatches found for ibin: {}".format(ibin))
            refined_extracted.append(ibin)
            continue     
        
        # If there are no methylation defined as at least one motif with a mean above 0.5: continue
        motif_features_for_bin = data.loc[data.index.isin(ibin), motifs]
        motif_features_for_bin = (motif_features_for_bin.mean() > 0.5).astype(int)
        if motif_features_for_bin.sum() == 0:
            logger.debug("No methylation pattern found for ibin: {}".format(ibin))
            refined_extracted.append(ibin)
            continue
        
        # get remaining eps values
        remaining_eps = [x for x in eps_values if x > float(ibin_eps)]
        
        # Create a copy of ibin to modify
        nbin = ibin.copy()
        tmp_bin = ibin.copy()
        
        for eps in remaining_eps:
            # Get cluster label for rbin (refined bin)
            rbin_label = results_df\
                .filter(pl.col("Contig").is_in(ibin))\
                .get_column(f"Cluster_Label_eps_{eps}")
            
            assert len(rbin_label.unique()) == 1, "rbin should have a unique cluster label"
            
            rbin_label = rbin_label[0]
            
            # Get rbin contigs
            rbin = unbinned_df\
                .filter(pl.col(f"Cluster_Label_eps_{eps}") == rbin_label)\
                .get_column("Contig")
            
            if len(rbin) == 0:
                continue
            
            # Add rbin to ibin
            tmp_bin.extend(rbin)
            # Remove duplicates
            tmp_bin = list(set(tmp_bin))
            
            rbin_pattern = methylation_comparison\
                .filter(pl.col("contig").is_in(tmp_bin))\
                .filter(pl.col("contig_compare").is_in(tmp_bin))
            
            # Check the rbin for mismatches
            if sum(rbin_pattern["sum_mismatches"] > 0):
                logger.debug("Mismatches found for eps: {}".format(eps))
                break
            
            nbin.extend(rbin)
            # Remove duplicates
            nbin = list(set(nbin))
        
        # Remove contigs from ubinned_df
        unbinned_df = unbinned_df\
            .filter(~pl.col("Contig").is_in(nbin))
        
        # Make sure all contigs are unique
        all_contigs = [item for sublist in refined_extracted for item in sublist]
        assert len(all_contigs) == len(set(all_contigs)), "Not all contigs are unique"
        refined_extracted.append(nbin)
        
        t_contigs_added += len(nbin) - len(ibin)
        
    all_contigs = [item for sublist in refined_extracted for item in sublist]
    assert len(all_contigs) == len(set(all_contigs)), "Not all contigs are unique"
    logger.info("Total contigs added: {}".format(t_contigs_added))
    return refined_extracted







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

    if len(features_data["motif"]) > 0:
        # Prepare DataFrame for saving results
        unbinned_df = pd.DataFrame({'Contig': contig_list})

        # Add cluster labels for each eps value to the DataFrame
        for eps_value in eps_values:
            unbinned_df[f'Cluster_Label_eps_{eps_value}'] = DBSCAN_results_dict[eps_value]

        unbinned_df = pl.DataFrame(unbinned_df)
        print(unbinned_df.shape)
        print("extracted length:", len([item for sublist in extracted for item in sublist]))
        extracted = expand_cluster_based_on_methylation(
            extracted,
            results_df,
            unbinned_df,
            eps_values,
            data,
            features_data,
            logger,
            args
        )

        print("refined extracted length:", len([item for sublist in extracted for item in sublist]))
        
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
