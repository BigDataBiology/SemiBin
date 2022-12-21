from SemiBin import naive_orffinder, fasta
import os
from os import path
from hypothesis import given, strategies as st

def test_find_orfs():
    orfs = naive_orffinder.find_orfs('ATG' + 'AAG' * 100 + 'TAA')
    assert orfs == [(0, 306, False)]

@given(seq=st.text(min_size=203, max_size=903, alphabet='ATGC'),
        pats=st.lists(min_size=1, max_size=7,
            elements=st.text(min_size=1, max_size=5, alphabet='ATGC')))
def test_findall(seq, pats):
    matches = naive_orffinder.findall(seq, pats)
    matches = set(matches)
    for i in range(len(seq)):
        is_match = any(seq[i:].startswith(p) for p in pats)
        assert (i in matches) == is_match


def test_run_naiveorf(tmpdir):
    ifile = 'test/single_sample_data/input.fasta.gz'

    os.makedirs(f'{tmpdir}/run1')
    oname1 = naive_orffinder.run_naiveorf(ifile,
                1,
                f'{tmpdir}/run1')

    os.makedirs(f'{tmpdir}/run2')
    oname2 = naive_orffinder.run_naiveorf(ifile,
                2,
                f'{tmpdir}/run2')
    seqs1 = list(fasta.fasta_iter(oname1))
    seqs2 = list(fasta.fasta_iter(oname2))
    assert seqs1 == seqs2

    ## Below are known values for the test data
    assert len(seqs1) == 3550

    seqs1.sort()
    assert seqs1[35] == ('g1k_0_40', 'LAICMFFTTNNIIKPEPNICPEIPILNNPALKAKINDNAAIPKIATFVK*')
    assert seqs1[125] == ('g1k_1_37', 'VENFLLNYQQIVVIIPLRMRQELLLKDVLEKLIPSQFHHILGLVRHQQLD*')

