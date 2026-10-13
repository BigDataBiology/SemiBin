from __future__ import annotations
import argparse
import logging
import os
from os import path
import math
import shutil
import tempfile
from typing import TYPE_CHECKING

from .utils import write_bins, log_binning_stats
from .markers import estimate_seeds

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd
    import torch

# This is the default in the igraph package
NR_INFOMAP_TRIALS = 10

def cal_kl(means: np.ndarray, variances: np.ndarray, use_ne: str = 'auto') -> np.ndarray:
    '''Pairwise KL divergence between contig coverage distributions

    Coverage of each contig is modelled as a normal distribution. Values are
    clipped to the [1e-6, 1-1e-6] range.

    Parameters
    ----------
    means : np.ndarray (1D)
        Mean coverage of each contig
    variances : np.ndarray (1D)
        Variance of the coverage of each contig
    use_ne : str, optional
        Whether to use numexpr (which is faster and uses less memory):
        'auto' (use it if available), 'yes' (fail if not available), or 'no'

    Returns
    -------
    kl : np.ndarray (2D)
        Matrix of size N x N, where N is the number of contigs
    '''
    # A naive implementation creates a lot of copies of what can
    # become large matrices
    import numpy as np

    means = np.clip(means, 1e-6, None)
    variances = np.clip(variances, 1.0, None)

    m1 = means.reshape(1, len(means))
    m2 = means.reshape(len(means), 1)

    v1 = variances.reshape(1, len(variances))
    v2 = variances.reshape(len(variances), 1)


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



def run_embed_infomap(logger: logging.Logger,
                      model: torch.nn.Module,
                      data: pd.DataFrame, * ,
                      device: torch.device,
                      max_edges: int,
                      max_node: float,
                      is_combined: bool,
                      n_sample: int,
                      contig_dict: dict[str, str],
                      num_process: int,
                      random_seed: int | None) -> tuple[np.ndarray, np.ndarray]:
    """
    Cluster contigs into bins

    Contigs are embedded with the model and a k-nearest-neighbour graph is
    built (keeping only edges that are also present in the k-NN graph of the
    input features). If coverage is not part of the model input, edge weights
    are additionally scaled by coverage similarity. The graph is then clustered
    with Infomap.

    Parameters
    ----------
    logger : logging.Logger
    model : torch.nn.Module
        Trained embedding model
    data : pd.DataFrame
        Contig features (indexed by contig name)
    device : torch.device
    max_edges : int
        max edges of one contig considered in binning
    max_node : float
        max fraction of contigs considered in binning
    is_combined : bool
        Whether coverage is part of the model input
    n_sample : int
        Number of samples
    contig_dict : dict[str, str]
        contig name -> sequence
    num_process : int
        Number of processes to use
    random_seed : int or None
        Currently unused

    Returns
    -------
    embedding : np.ndarray
        Embedding of each contig
    contig_labels : np.ndarray
        Cluster label of each contig
    """
    from igraph import Graph
    from sklearn.neighbors import kneighbors_graph
    from scipy import sparse
    import torch
    import numpy as np
    from .utils import norm_abundance

    train_data_input = data.values[:, 0:136] if not is_combined else data.values

    if is_combined:
        if norm_abundance(train_data_input):
            from sklearn.preprocessing import normalize
            norm = np.sum(train_data_input, axis=0)
            train_data_input = train_data_input / norm
            train_data_input = normalize(train_data_input, axis=1, norm='l1')

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

    edges = [(x,y) for x,y in zip(X, Y, strict=True)]
    edge_weights = V[above_diag]

    logger.debug(f'Number of edges in clustering graph: {len(edges)}')

    g = Graph()
    g.add_vertices(np.arange(num_contigs))
    g.add_edges(edges)
    length_weight = np.array([len(contig_dict[name]) for name in data.index])
    logger.debug(f'Running infomap with {num_process} processes...')
    result = g.community_infomap(
                edge_weights=edge_weights,
                vertex_weights=length_weight,
                trials=NR_INFOMAP_TRIALS)
    contig_labels = np.zeros(shape=num_contigs, dtype=int)

    for cluster_ix, members in enumerate(result):
        contig_labels[members] = cluster_ix
    return embedding, contig_labels


def recluster_bins(logger: logging.Logger,
                   data: pd.DataFrame, *,
                   n_sample: int,
                   embedding: np.ndarray,
                   is_combined: bool,
                   contig_labels: np.ndarray,
                   minfasta: int,
                   contig_dict: dict[str, str],
                   binned_length: int,
                   orf_finder: str,
                   num_process: int,
                   random_seed: int | None) -> np.ndarray:
    '''Split bins that appear to contain multiple genomes

    Single-copy marker genes are used to estimate how many genomes each bin
    contains. Bins with more than one estimated genome are split with
    (length-weighted) k-means, seeded with the contigs containing the marker
    genes.

    Parameters
    ----------
    logger : logging.Logger
    data : pd.DataFrame
        Contig features (indexed by contig name)
    n_sample : int
        Number of samples
    embedding : np.ndarray
        Embedding of each contig (as returned by ``run_embed_infomap``)
    is_combined : bool
        Whether coverage is part of the model input. If False, the mean
        coverage is appended to the embedding before k-means
    contig_labels : np.ndarray
        Initial bin label of each contig
    minfasta : int
        Minimum bin size (in basepairs); smaller bins are not reclustered
    contig_dict : dict[str, str]
        contig name -> sequence
    binned_length : int
        Minimum contig length
    orf_finder : str
        ORF finder to use
    num_process : int
        Number of processes to use
    random_seed : int or None
        Random seed for k-means

    Returns
    -------
    contig_labels_reclustered : np.ndarray
        New bin label of each contig (all labels are non-negative)
    '''
    from sklearn.cluster import KMeans
    import numpy as np
    from collections import defaultdict
    if not is_combined:
        depth = data.values[:, 136:len(data.values[0])].astype(np.float32)
        mean_index = [2 * sample_ix for sample_ix in range(n_sample)]
        depth_mean = depth[:, mean_index]
        abun_scale = np.ceil(depth_mean.mean(axis = 0) / 100) * 100
        depth_mean = depth_mean / abun_scale
        scaling = np.mean(np.abs(embedding)) / np.mean(depth_mean)
        base = 10
        weight = 2 * base * math.ceil(scaling / base)
        embedding_with_depth = np.concatenate(
            (embedding, depth_mean * weight), axis=1)
    else:
        embedding_with_depth = embedding


    total_size = defaultdict(int)
    for name, label in zip(data.index, contig_labels, strict=True):
        total_size[label] += len(contig_dict[name])
    with tempfile.TemporaryDirectory() as tdir:
        cfasta = os.path.join(tdir, 'concatenated.fna')
        with open(cfasta, 'wt') as concat_out:
            for name, bin_ix in zip(data.index, contig_labels, strict=True):
                if total_size[bin_ix] < minfasta:
                    continue
                concat_out.write(f'>bin{bin_ix:06}.{name}\n')
                concat_out.write(contig_dict[name] + '\n')

        seeds = estimate_seeds(
            cfasta,
            binned_length,
            num_process,
            multi_mode=True,
            orf_finder=orf_finder)
            # we cannot bypass the orf_finder here, because of the renaming of the contigs
        if not seeds:
            logger.warning('No bins found in the concatenated fasta file.')
            return contig_labels

    name2ix = {name:ix for ix,name in enumerate(data.index)}
    contig_labels_reclustered = np.empty_like(contig_labels)
    contig_labels_reclustered.fill(-1)
    next_label = 0
    for bin_ix in range(contig_labels.max() + 1):
        bin_seeds = seeds.get(f'bin{bin_ix:06}', [])
        n_subbins = len(bin_seeds)

        if n_subbins > 1 and total_size[bin_ix] >= minfasta:
            contig_indices = [i for i,ell in enumerate(contig_labels) if ell == bin_ix]
            bin_features = embedding_with_depth[contig_indices]

            seed_index = [name2ix[s] for s in bin_seeds]
            length_weight = np.array(
                [len(contig_dict[name])
                    for name,ell in zip(data.index, contig_labels, strict=True)
                        if ell == bin_ix])
            seeds_embedding = embedding_with_depth[seed_index]
            kmeans = KMeans(
                n_clusters=n_subbins,
                init=seeds_embedding,
                n_init=1,
                random_state=random_seed)
            kmeans.fit(bin_features, sample_weight=length_weight)
            for contig_ix, label in zip(contig_indices, kmeans.labels_, strict=True):
                contig_labels_reclustered[contig_ix] = next_label + label
            next_label += n_subbins
        else:
            contig_labels_reclustered[contig_labels == bin_ix] = next_label
            next_label += 1
    assert contig_labels_reclustered.min() >= 0
    return contig_labels_reclustered


def cluster(logger: logging.Logger,
            model: torch.nn.Module,
            data: pd.DataFrame,
            device: torch.device,
            is_combined: bool,
            n_sample: int,
            out: str,
            contig_dict: dict[str, str], *,
            args: argparse.Namespace,
            binned_length: int,
            minfasta: int) -> None:
    """
    Cluster contigs into bins (short-read mode) and write the output

    Runs ``run_embed_infomap`` and then (unless disabled with
    ``--no-recluster``) ``recluster_bins``. Bins are written to subdirectories
    of ``out`` along with ``contig_bins.tsv`` and a TSV file with information
    on the bins.

    Parameters
    ----------
    logger : logging.Logger
    model : torch.nn.Module
        Trained embedding model
    data : pd.DataFrame
        Contig features (indexed by contig name)
    device : torch.device
    is_combined : bool
        Whether coverage is part of the model input
    n_sample : int
        Number of samples
    out : str
        Output directory
    contig_dict : dict[str, str]
        contig name -> sequence
    args : argparse.Namespace
        Command line arguments (for ``max_edges``, ``max_node``,
        ``recluster``, ``write_pre_reclustering_bins``, output options...)
    binned_length : int
        Minimum contig length
    minfasta : int
        Minimum bin size (in basepairs) for a bin to be output
    """
    import pandas as pd
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
            log_binning_stats(logger, bin_files, contig_dict, binned_length)
        n_pre_bins = len(bin_files)
    else:
        from collections import defaultdict
        total_size = defaultdict(int)
        for name, label in zip(data.index, contig_labels, strict=True):
            total_size[label] += len(contig_dict[name])
        n_pre_bins = sum((total_size[bin_ix] >= minfasta for bin_ix in range(contig_labels.max() + 1)))


    if args.recluster:
        logger.info(f'Number of bins prior to reclustering: {n_pre_bins}')
        logger.debug('Reclustering...')

        if n_pre_bins == 0:
            import numpy as np
            logger.warning('No bins were created. Please check your input data.')
            contig_labels_reclustered = np.full(len(data.index), fill_value=-1, dtype=int)
        else:
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
        log_binning_stats(logger, outputs, contig_dict, binned_length)
    logger.info('Binning finished')

