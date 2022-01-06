from SemiBin.generate_coverage import calculate_coverage
import pandas as pd
from pandas.testing import assert_frame_equal


def test_cov():

    test_data = [
        'k141_63080\t0\t6\t1\n',
        'k141_63080\t6\t30\t2\n',
        'k141_63080\t30\t138\t3\n',
        'k141_63080\t138\t144\t2\n',
        'k141_63080\t144\t174\t1\n',
        'k141_63080\t174\t177\t2\n',
        'k141_63080\t177\t257\t3\n',
        'k141_63080\t257\t304\t2\n',
        'k141_63080\t304\t315\t1\n',
        'k141_0\t0\t45\t1\n',
        'k141_0\t45\t88\t2\n',
        'k141_0\t88\t138\t3\n',
        'k141_0\t138\t159\t2\n',
        'k141_0\t159\t183\t3\n',
        'k141_0\t183\t225\t2\n',
        'k141_0\t225\t227\t1\n',
        'k141_0\t227\t297\t3\n',
        'k141_0\t297\t338\t2\n',
        'k141_0\t338\t365\t1\n',
    ]

    cov, cov_split = calculate_coverage(
        test_data,
        'test_data',
        must_link_threshold=0,
        is_combined=True,
        contig_threshold=0)
    cov.columns = ['cov']
    cov_split.columns = ['cov']
    assert_frame_equal(cov, pd.DataFrame([2.581818, 2.627907], index=[
                       'k141_63080', 'k141_0'], columns=['cov']))
    assert_frame_equal(cov_split, pd.DataFrame([2.609756, 2.554217, 2.682243, 2.574074], index=[
                       'k141_63080_1', 'k141_63080_2', 'k141_0_1', 'k141_0_2'], columns=['cov']))

    cov = calculate_coverage(
        test_data,
        'test_data',
        must_link_threshold=0,
        is_combined=False,
        contig_threshold=0)
    cov.columns = ['mean', 'var']
    assert_frame_equal(cov, pd.DataFrame([[2.581818, 0.606942], [
                       2.627907, 0.252244]], index=['k141_63080', 'k141_0'], columns=['mean', 'var']))
