"""
https://github.com/BinPro/CONCOCT/blob/develop/scripts/fasta_to_features.py
"""
from itertools import product
from Bio import SeqIO
from itertools import tee
from collections import Counter, OrderedDict
from Bio.SeqRecord import SeqRecord


def window(seq, n):
    els = tee(seq, n)
    for i, el in enumerate(els):
        for _ in range(i):
            next(el, None)
    return zip(*els)


def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC", repeat=kmer_len):
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash, counter


def generate_kmer_features_from_fasta(
        fasta_file, length_threshold, kmer_len, split=False, threshold=0):
    import numpy as np
    import pandas as pd
    kmer_dict, nr_features = generate_feature_mapping(kmer_len)
    composition_d = OrderedDict()
    contig_lengths = OrderedDict()
    if not split:
        seq_list = list(SeqIO.parse(fasta_file, "fasta"))
    else:
        seq_list = []
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if len(seq_record) >= threshold:
                half = int(len(seq_record.seq) / 2)
                rec1 = SeqRecord(
                    seq_record.seq[0:half], id=seq_record.id + '_1', description='')
                seq_list.append(rec1)
                rec2 = SeqRecord(seq_record.seq[half:len(
                    seq_record.seq)], id=seq_record.id + '_2', description='')
                seq_list.append(rec2)
    for seq in seq_list:
        seq_len = len(seq)
        if seq_len <= length_threshold:
            continue
        contig_lengths[seq.id] = seq_len
        kmers = [
            kmer_dict[kmer_tuple]
            for kmer_tuple
            in window(str(seq.seq).upper(), kmer_len)
            if kmer_tuple in kmer_dict
        ]
        kmers.append(nr_features - 1)
        composition_v = np.bincount(np.array(kmers, dtype=np.int64))
        composition_v[-1] -= 1
        composition_d[seq.id] = composition_v
    df = pd.DataFrame.from_dict(composition_d, orient='index', dtype=float)

    df = df.apply(lambda x: x + 1e-5)
    df = df.div(df.sum(axis=1), axis=0)
    return df
