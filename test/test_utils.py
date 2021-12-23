from SemiBin.utils import get_must_link_threshold, get_marker
from hypothesis import given, strategies as st

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


def test_get_marker():
    rs = get_marker('test/data/bin.230.fa.hmmout')
    assert rs == ['k119_487808', 'k119_606328']
    rs = get_marker('test/data/bin.96.fa.hmmout')
    assert rs == ['k119_268294', 'k119_337646']

