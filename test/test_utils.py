from SemiBin.utils import get_must_link_threshold, split_data, n50_l50, maybe_crams2bams, norm_abundance, normalize_kmer_motif_features
from SemiBin import utils
from hypothesis import given, strategies as st
import numpy as np
import pandas as pd


def slow_get_must_link_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    basepair_sum = 0
    threshold = 0
    whole_len = np.sum(contig_len)
    contig_len.sort(reverse=True)
    index = 0
    while(basepair_sum / whole_len < 0.98):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    return np.clip(threshold, 4000, None)


@given(lens=st.lists(st.integers(min_value=1, max_value=(1<<32)), min_size=1))
def test_get_must_link_threshold(lens):
    assert get_must_link_threshold(lens) == slow_get_must_link_threshold(lens)


def test_split_data():
    '''Regression test for #68

    https://github.com/BigDataBiology/SemiBin/issues/68
    '''
    import pandas as pd
    data = pd.DataFrame([[3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02],
           [3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02],
           [3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02],
           [3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02, 3.02367003e-02, 3.02367003e-02,
            3.02367003e-02, 3.02367003e-02],
           [1.00000000e-05, 1.00000000e-05, 1.00000000e-05, 1.00000000e-05,
            1.00000000e-05, 1.00000000e-05, 1.00000000e-05, 1.00000000e-05,
            1.00000000e-05, 1.00000000e-05]],
           index=[
               'S1_2341+340_METAG:g1k_5_1',
               'S1_2341+340_METAG:g2k_4_1',
               'S2_2341+340_METAG:g1k_2_1',
               'S2_2341+340_METAG:g1k_3_1',
               'S2_2341+340_METAG:g1k_5_1'],
            columns=[
               'sorted10.sam_cov',
               'sorted1.sam_cov',
               'sorted2.sam_cov',
               'sorted3.sam_cov',
               'sorted4.sam_cov',
               'sorted5.sam_cov',
               'sorted6.sam_cov',
               'sorted7.sam_cov',
               'sorted8.sam_cov',
               'sorted9.sam_cov'],
            ).reset_index().rename(columns={'index':'contig_name'})

    split1 = split_data(data, 'S1_2341+340_METAG', ':')
    assert len(split1) == 2
    split2 = split_data(data, 'S2_2341+340_METAG', ':')
    assert len(split2) == 3


@given(sizes=st.lists(st.integers(min_value=2, max_value=(1<<32)), min_size=1))
def test_n50_l50(sizes):
    n50,l50 = n50_l50(sizes)
    sizes = np.array(sizes)
    sizes.sort()
    sizes = sizes[::-1]
    assert l50 <= len(sizes)
    assert np.sum(sizes[:l50]) >= np.sum(sizes)//2
    assert np.sum(sizes[:l50-1]) < np.sum(sizes)//2
    assert n50 == sizes[l50-1]

def test_extract_bams(tmpdir):
    single_sample_input = 'test/single_sample_data'
    rs = maybe_crams2bams([f'{single_sample_input}/input.cram'],
             f'{single_sample_input}/input.fasta',
             2,
             tmpdir)
    assert len(rs) == 1
    assert rs[0].endswith('.bam')
    assert rs[0].startswith(str(tmpdir))

def test_normalize_kmer_motif_features():

    # Test case 1: Normalizing a 3x3 array
    data = {
        "kmer1": [1, 2, 3],
        "kmer2": [4, 5, 6],
        "kmer3": [7, 8, 9]
    }

    df = pd.DataFrame(data)
    df = df.values
    df_norm, _ = normalize_kmer_motif_features(df, df)

    
    # Expected normalized and concatenated DataFrame
    expected_df_norm = np.array([
        [0., 0., 0.],
        [0.5, 0.5, 0.5],
        [1., 1., 1.]
    ])
    assert np.allclose(df_norm, expected_df_norm), f"Expected {expected_df_norm}, but got {df_norm}"

    # Test case 2: Normalizing train_data with a different train_data_split
    df_split = pd.DataFrame({
        "kmer1": [1, 2],
        "kmer2": [4, 5],
        "kmer3": [7, 8]
    })
    df_split = df_split.values

    df_norm, df_split_norm = normalize_kmer_motif_features(df, df_split)
    
    # We normalize the split df values based on the full df. That is why we expect 0.5 
    expected_df_split_norm = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5]
    ])
    assert np.allclose(df_norm, expected_df_norm), f"Expected {expected_df_norm}, but got {df_norm}"
    assert np.allclose(df_split_norm, expected_df_split_norm), f"Expected {expected_df_split_norm}, but got {df_split_norm}"

def test_norm_abundance():
    features = {}
    features["kmer"] = [str(i) for i in range(136)]
    features["motif"] = []
    features["motif_present"] = []
    features["depth"] = []
    
    kmer_cols = [str(i) for i in range(136)]

    df = pd.DataFrame(np.random.randn(10, 136), columns=kmer_cols)
    assert not norm_abundance(df, features)

    features["depth"] = ["depth_1"]
    df = pd.DataFrame(np.random.randn(12, 137), columns=kmer_cols + features["depth"])
    assert not norm_abundance(df, features)

    features["depth"] = ["depth_1", "depth_2", "depth_3", "depth_4"]
    df = pd.DataFrame(np.random.randn(12, 140), columns=kmer_cols + features["depth"])
    assert not norm_abundance(df, features)

    features["depth"] = ["depth_1", "depth_2"]
    df = pd.DataFrame(np.abs(np.random.randn(12, 138))*2, columns=kmer_cols + features["depth"])
    assert not norm_abundance(df, features)

    features["depth"] = ["depth_1", "depth_2"]
    df = pd.DataFrame(np.abs(np.random.randn(12, 138))*4, columns=kmer_cols + features["depth"])
    assert not norm_abundance(df, features)

    features["depth"] = [f"depth_{i}" for i in range(12)]
    df = pd.DataFrame(np.abs(np.random.randn(12, 148)), columns=kmer_cols + features["depth"])
    assert norm_abundance(df, features)

    features["depth"] = [f"depth_{i}" for i in range(20)]
    df = pd.DataFrame(       np.random.randn(12, 156) , columns=kmer_cols + features["depth"])
    assert norm_abundance(df, features)

    features["depth"] = [f"depth_{i}" for i in range(28)]
    df = pd.DataFrame(       np.random.randn(12, 164) , columns=kmer_cols + features["depth"])
    assert norm_abundance(df, features)

    df = pd.DataFrame(np.abs(np.random.randn(12, 164)), columns=kmer_cols + features["depth"])
    assert norm_abundance(df, features)


def test_load_fasta():
    c_min_len, ml_thresh, inputs = utils.load_fasta('test/train_data/input.fasta', .01)
    assert c_min_len == 1000
