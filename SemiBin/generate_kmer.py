# Adapted from https://github.com/BinPro/CONCOCT/blob/develop/scripts/fasta_to_features.py
from itertools import product
from collections import OrderedDict

from .fasta import fasta_iter

def generate_feature_mapping(kmer_len):
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
        fasta_file, length_threshold, kmer_len, split=False, split_threshold=0):
    import numpy as np
    import pandas as pd
    def seq_list():
        for h, seq in fasta_iter(fasta_file):
            if not split:
                yield h, seq
            elif len(seq) >= split_threshold:
                half = len(seq) // 2
                yield (h + '_1', seq[:half])
                yield (h + '_2', seq[half:])

    kmer_dict, nr_features = generate_feature_mapping(kmer_len)
    composition = OrderedDict()
    for h, seq in seq_list():
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
