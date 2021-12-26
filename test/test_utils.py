from SemiBin.utils import get_must_link_threshold, get_marker
from hypothesis import given, strategies as st
from io import StringIO

def slow_get_must_link_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    import numpy as np
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


HMMOUT_TEST_SINGLE = '''#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
k119_224042_6410_6778_+ -            122 TIGR00855            TIGR00855    125     3e-39  129.7   6.7   1   1   4.4e-42   3.3e-39  129.6   6.7     2   125     2   122     1   122 0.92 -
#
# Program:         hmmsearch
# Version:         3.3.2 (Nov 2020)
# Pipeline mode:   SEARCH
# Query file:      /home/luispedro/.miniconda3/envs/py3.9/lib/python3.9/site-packages/SemiBin-0.4.0-py3.9.egg/SemiBin/marker.hmm
# Target file:     restart/bin.50.fa.frag.faa
# Option settings: hmmsearch --domtblout restart/bin.50.fa.hmmout --cut_tc --cpu 4 /home/luispedro/.miniconda3/envs/py3.9/lib/python3.9/site-packages/SemiBin-0.4.0-py3.9.egg/SemiBin/marker.hmm restart/bin.50.fa.frag.faa
# Current dir:     /home/luispedro/Sync/work/SemiBin.rec
# [ok]
'''

def test_get_marker():
    rs = get_marker('test/data/bin.230.fa.hmmout')
    assert rs == ['k119_487808', 'k119_606328']
    rs = get_marker('test/data/bin.96.fa.hmmout')
    assert rs == ['k119_268294', 'k119_337646']

    # This caused a crash in a previous version of SemiBin
    rs = get_marker(StringIO(HMMOUT_TEST_SINGLE))
    assert rs == ['k119_224042']

