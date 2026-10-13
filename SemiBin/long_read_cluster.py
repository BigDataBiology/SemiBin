from __future__ import annotations
import argparse
import logging
import os
import shutil
from typing import TYPE_CHECKING
import torch
import numpy as np
from .utils import write_bins, log_binning_stats
from .markers import estimate_seeds, get_marker, N_SINGLE_COPY_MARKERS
from sklearn.cluster import dbscan
from sklearn.neighbors import kneighbors_graph
from collections import defaultdict

if TYPE_CHECKING:
    import pandas as pd


def get_best_bin(results: list[list[int]],
                 contig_to_marker: dict[str, list[str]],
                 namelist: list[str],
                 contig_dict: dict[str, str],
                 minfasta: int) -> list[str] | None:
    '''Select the best bin across several clustering results

    Bins are scored by an F1-like combination of completeness and
    contamination, estimated from single-copy marker genes. The contamination
    limit is progressively relaxed (from 10% up to 100%) until some bin
    qualifies. Ties are broken in favour of the smaller bin.

    Parameters
    ----------
    results : list of list of int
        Each element is a clustering: one label per contig in ``namelist``
        (-1 means unclustered)
    contig_to_marker : dict[str, list[str]]
        contig name -> marker genes found in it
    namelist : list of str
        Contig names
    contig_dict : dict[str, str]
        contig name -> sequence
    minfasta : int
        Minimum bin size (in basepairs)

    Returns
    -------
    best_bin : list of str or None
        Contigs in the best bin, or None if no bin of at least ``minfasta``
        basepairs contains any marker genes
    '''

    # There is room for improving the loop below to avoid repeated computation
    # but it runs very fast in any case
    for max_contamination in [0.1, 0.2, 0.3, 0.4, 0.5, 1]:
        max_F1 = 0
        weight_of_max = 1e9
        max_bin = None

        for labels in results:
            bin_contigs = defaultdict(list)
            for label, name in zip(labels, namelist, strict=True):
                if label != -1:
                    bin_contigs[label].append(name)
            for bin_contig in bin_contigs.values():
                cur_weight = sum(len(contig_dict[contig]) for contig in bin_contig)
                if cur_weight < minfasta:
                    continue
                marker_list = []
                for contig in bin_contig:
                    marker_list.extend(contig_to_marker[contig])
                if len(marker_list) == 0:
                    continue
                recall = len(set(marker_list)) / N_SINGLE_COPY_MARKERS
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

def cluster_long_read(logger: logging.Logger,
                      model: torch.nn.Module,
                      data: pd.DataFrame,
                      device: torch.device,
                      is_combined: bool,
                      n_sample: int,
                      out: str,
                      contig_dict: dict[str, str], *,
                      binned_length: int,
                      args: argparse.Namespace,
                      minfasta: int) -> None:
    '''Cluster contigs into bins (long-read mode) and write the output

    Contigs are embedded with the model and clustered with DBSCAN using a
    range of ``eps`` values. Bins are then greedily extracted from this
    ensemble of clusterings: at each step, the best bin (according to
    single-copy marker genes, see ``get_best_bin``) is selected and its
    contigs are removed from consideration.

    Bins are written to ``<out>/output_bins`` along with ``bins_info.tsv`` and
    ``contig_bins.tsv``.

    Parameters
    ----------
    logger : logging.Logger
    model : torch.nn.Module
        Trained embedding model
    data : pd.DataFrame
        Contig features (indexed by contig name)
    device : torch.device
    is_combined : bool
        Whether coverage is part of the model input. If False, the log of the
        mean coverage is appended to the embedding
    n_sample : int
        Number of samples
    out : str
        Output directory
    contig_dict : dict[str, str]
        contig name -> sequence
    binned_length : int
        Minimum contig length
    args : argparse.Namespace
        Command line arguments (for ``num_process``, ``orf_finder``, output
        options...)
    minfasta : int
        Minimum bin size (in basepairs) for a bin to be output
    '''
    import pandas as pd
    from .utils import norm_abundance
    contig_list = data.index.tolist()
    if not is_combined:
        train_data_input = data.values[:, 0:136]
    else:
        train_data_input = data.values
        if norm_abundance(train_data_input):
            from sklearn.preprocessing import normalize
            norm = np.sum(train_data_input, axis=0)
            train_data_input = train_data_input / norm
            train_data_input = normalize(train_data_input, axis=1, norm='l1')

    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()

    length_weight = np.array(
        [len(contig_dict[name]) for name in contig_list])

    if not is_combined:
        depth = data.values[:, 136:len(data.values[0])].astype(np.float32)
        mean_index = [2 * sample_ix for sample_ix in range(n_sample)]
        depth = depth[:, mean_index]
        embedding_with_depth = np.concatenate((embedding, np.log(np.clip(depth, 1e-6, None))), axis=1)
    else:
        embedding_with_depth = embedding

    import tempfile
    with tempfile.TemporaryDirectory() as tdir:
        cfasta = os.path.join(tdir, 'concatenated.fna')
        with open(cfasta, 'wt') as concat_out:
            for h in contig_list:
                concat_out.write(f'>{h}\n{contig_dict[h]}\n')
            #  This is needed to create the markers.hmmout file
            estimate_seeds(
                cfasta,
                binned_length,
                args.num_process,
                output=out,
                orf_finder=args.orf_finder,
                prodigal_output_faa=args.prodigal_output_faa)

            contig2marker = get_marker(os.path.join(out, 'markers.hmmout'), orf_finder=args.orf_finder,
                                       min_contig_len=binned_length, fasta_path=cfasta, contig_to_marker=True)

    output_bin_path = os.path.join(out, 'output_bins')
    if os.path.exists(output_bin_path):
        logger.warning(f'Previous output directory `{output_bin_path}` found. Over-writing it.')
        shutil.rmtree(output_bin_path)
    os.makedirs(output_bin_path, exist_ok=True)

    logger.debug('Running DBSCAN.')
    dist_matrix = kneighbors_graph(
        embedding_with_depth,
        n_neighbors=min(200, embedding_with_depth.shape[0] - 1),
        mode='distance',
        p=2,
        n_jobs=args.num_process)

    dbscan_results = []
    for eps_value in [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55]:
        _, labels = dbscan(dist_matrix,
                eps=eps_value, min_samples=5, n_jobs=args.num_process, metric='precomputed',
                sample_weight=length_weight)
        dbscan_results.append(labels.tolist())

    logger.debug('Integrating results.')

    extracted = []
    while sum(len(contig_dict[contig]) for contig in contig_list) >= minfasta:
        if len(contig_list) == 1:
            extracted.append(contig_list)
            break

        max_bin = get_best_bin(dbscan_results,
                                contig2marker,
                                contig_list,
                                contig_dict,
                                minfasta)
        if not max_bin:
            break

        extracted.append(max_bin)
        for clustered in max_bin:
            ix = contig_list.index(clustered)
            contig_list.pop(ix)
            for r in dbscan_results:
                r.pop(ix)

    contig2bin = {}
    for bin_ix, bin_contigs in enumerate(extracted):
        for contig in bin_contigs:
            contig2bin[contig] = bin_ix
    namelist = data.index.tolist()
    contig_labels = [contig2bin.get(c, -1) for c in namelist]
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
    log_binning_stats(logger, written, contig_dict, binned_length)
    logger.info('Finished binning.')
