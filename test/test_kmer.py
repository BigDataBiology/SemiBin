from SemiBin.generate_kmer import generate_feature_mapping, generate_kmer_features_from_fasta
from os import path
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import numpy as np

def test_generate_feature_mapping():
    fm,n = generate_feature_mapping(4)
    assert n == len(set(fm.values()))


def test_kmer(tmpdir):
    fasta_file = path.join(tmpdir, 'contig.fna')

    with open(fasta_file, 'wt') as out_fasta:
        out_fasta.write(f'>example\nGCATGAGACAGTGTTTGCTGCAATACTCCCTGAGGCGCTTGCATATTATGAGAAGAATGCTCCCAAGGGAGAGTGCGTGATCGTTATAGAGGGAAAGAGCAGGCTTGAGATTCGGGAAGAAGAGAAAGCACAGTGGGAACAGCTTACCATTGAAGAACATATGGAACGCTATCTCTCCGGCGGTATGGAGAAAAAGGAAGCTATGAAGCGGGTGGCAAAGGATCGGGGCGTAAGTAAAAGAGATATTTATCAAGCGCTTCTTTAAACCTTCATATGAATGAAATAAAAACAGGTATCGGGTGTGTTC\n')

    kmer = generate_kmer_features_from_fasta(
        fasta_file, 0, kmer_len=1, split=False, split_threshold=0)
    assert_frame_equal(kmer, pd.DataFrame(
        np.array([0.547231, 0.452769]).reshape(1, 2), index=['example']))

    kmer_split = generate_kmer_features_from_fasta(
        fasta_file, 0, kmer_len=1, split=True, split_threshold=0)
    assert_frame_equal(kmer_split, pd.DataFrame(
        [[0.516340, 0.483660], [0.577922, 0.422078]], index=['example_1', 'example_2']))

    kmer = generate_kmer_features_from_fasta(
        fasta_file, 0, kmer_len=2, split=False, split_threshold=0)
    assert_frame_equal(kmer,
                       pd.DataFrame(np.array([0.169935,
                                              0.078431,
                                              0.147059,
                                              0.07516340680934419,
                                              0.058823542868123246,
                                              0.117647,
                                              0.143791,
                                              0.104575,
                                              0.06535948844461162,
                                              0.03921570613865813]).reshape(1, 10),
                                    index=['example']))

    kmer_split = generate_kmer_features_from_fasta(
        fasta_file, 0, kmer_len=2, split=True, split_threshold=0)
    assert_frame_equal(kmer_split,
                       pd.DataFrame([[0.131579, 0.06578949619112093, 0.177632, 0.07236843923128998, 0.03947372403044471, 0.151316, 0.157895, 0.098684, 0.07894738227145903, 0.02631583795010661],
                                    [0.202614, 0.0915032735272722, 0.117647, 0.07843138664615251, 0.07843138664615251, 0.084967, 0.130719, 0.111111, 0.05228761288391314, 0.05228761288391314]],
                                                index=['example_1', 'example_2']))

def test_kmer_with_ns(tmpdir):
    fasta_file = path.join(tmpdir, 'contig.fna')

    with open(fasta_file, 'wt') as out_fasta:
        out_fasta.write(f'>example\n')
        out_fasta.write('AAAAAAAAAAAAAAANAAAAAAAAAAAA\n')

    kmer = generate_kmer_features_from_fasta(
        fasta_file, 0, kmer_len=4, split=False, split_threshold=0)
    kmer = kmer.squeeze()
    assert len(kmer) == 136
    assert kmer.max() > 0.95

