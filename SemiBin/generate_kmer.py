# Adapted from https://github.com/BinPro/CONCOCT/blob/develop/scripts/fasta_to_features.py
from __future__ import annotations
from itertools import product
from collections import OrderedDict
from collections.abc import Iterator
from typing import TYPE_CHECKING

from .fasta import fasta_iter

if TYPE_CHECKING:
    import pandas as pd

def generate_feature_mapping(kmer_len: int) -> tuple[dict[str, int], int]:
    '''Map every k-mer to a feature index, merging reverse complements

    A k-mer and its reverse complement are assigned the same index (canonical
    k-mers), so that features do not depend on which strand was assembled.

    Parameters
    ----------
    kmer_len : int
        Length of the k-mers

    Returns
    -------
    kmer_hash : dict[str, int]
        k-mer -> feature index
    nr_features : int
        Number of distinct features (i.e., canonical k-mers)
    '''
    BASE_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC", repeat=kmer_len):
        kmer = ''.join(kmer)
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[''.join(rev_compl)] = counter
            counter += 1
    return kmer_hash, counter


def generate_kmer_features_from_fasta(
        fasta_file: str,
        length_threshold: int,
        kmer_len: int,
        split: bool = False,
        split_threshold: int = 0) -> pd.DataFrame:
    '''Compute normalized k-mer frequency features for each contig

    Parameters
    ----------
    fasta_file : str
        Input FASTA file
    length_threshold : int
        Sequences shorter than this (after splitting, if ``split`` is used) are
        skipped
    kmer_len : int
        Length of the k-mers (SemiBin uses 4)
    split : bool, optional
        If True, each contig of at least ``split_threshold`` bp is split into
        two halves (named ``<contig>_1`` and ``<contig>_2``), which are used to
        generate must-link constraints. Shorter contigs are dropped.
    split_threshold : int, optional
        Minimum contig length to be split (only used if ``split`` is True)

    Returns
    -------
    features : pd.DataFrame
        One row per (possibly split) contig, one column per canonical k-mer.
        Each row is a frequency vector (sums to 1), with a small pseudocount
        added.
    '''
    import numpy as np
    import pandas as pd
    def iter_sequences() -> Iterator[tuple[str, str]]:
        for h, seq in fasta_iter(fasta_file):
            if not split:
                yield h, seq
            elif len(seq) >= split_threshold:
                half = len(seq) // 2
                yield (h + '_1', seq[:half])
                yield (h + '_2', seq[half:])

    kmer_dict, nr_features = generate_feature_mapping(kmer_len)
    composition = OrderedDict()
    for h, seq in iter_sequences():
        if len(seq) < length_threshold:
            continue
        norm_seq = str(seq).upper()
        kmers = [kmer_dict[norm_seq[i:i+kmer_len]]
                for i in range(len(norm_seq) - kmer_len + 1)
                if norm_seq[i:i+kmer_len] in kmer_dict] # ignore kmers with non-canonical bases
        composition[h] = np.bincount(np.array(kmers, dtype=np.int64), minlength=nr_features)
    df = pd.DataFrame.from_dict(composition, orient='index', dtype=float)

    df = df.apply(lambda x: x + 1e-5)
    df = df.div(df.sum(axis=1), axis=0)
    return df
