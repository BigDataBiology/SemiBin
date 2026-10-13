import json
import os
import shutil
from unittest import mock

from SemiBin import markers
from SemiBin.markers import get_marker, estimate_seeds
from io import StringIO


HMMOUT_TEST_SINGLE = '''#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
k119_224042_1 -            122 TIGR00855            TIGR00855    125     3e-39  129.7   6.7   1   1   4.4e-42   3.3e-39  129.6   6.7     2   125     2   122     1   122 0.92 -
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
    rs = get_marker('test/data/bin.230.fa.hmmout', orf_finder='fraggenescan')
    assert rs == ['k119_487808', 'k119_606328']
    rs = get_marker('test/data/bin.96.fa.hmmout', orf_finder='fraggenescan')
    assert rs == ['k119_268294', 'k119_337646']

    # This caused a crash in a previous version of SemiBin
    rs = get_marker(StringIO(HMMOUT_TEST_SINGLE), orf_finder='prodigal')
    assert rs == ['k119_224042']

# Same hit as above, but only covering hmm positions 2-10 (of 125), so it is
# removed by the alignment-coverage filter
HMMOUT_TEST_LOW_COVERAGE = HMMOUT_TEST_SINGLE.replace(
        '    2   125     2   122     1   122',
        '    2    10     2   122     1   122')

def test_get_marker_all_filtered():
    rs = get_marker(StringIO(HMMOUT_TEST_LOW_COVERAGE), orf_finder='prodigal')
    assert rs == []

    rs = get_marker(StringIO(HMMOUT_TEST_LOW_COVERAGE.replace('k119_224042_1', 'bin000001.k119_224042_1')),
                    orf_finder='prodigal', multi_mode=True)
    assert rs == {}

    rs = get_marker(StringIO(HMMOUT_TEST_LOW_COVERAGE), orf_finder='prodigal', contig_to_marker=True)
    assert len(rs) == 0
    assert rs['k119_224042'] == []


def test_get_marker_all_filtered_by_length(tmp_path):
    fasta = tmp_path / 'contigs.fa'
    fasta.write_text('>k119_224042\n' + 'A' * 500 + '\n')
    rs = get_marker(StringIO(HMMOUT_TEST_SINGLE), fasta_path=str(fasta),
                    min_contig_len=1000, orf_finder='prodigal')
    assert rs == []


def test_get_marker_no_hits():
    no_hits = '\n'.join(ell for ell in HMMOUT_TEST_SINGLE.splitlines() if ell.startswith('#'))
    assert get_marker(StringIO(no_hits), orf_finder='prodigal') == []
    assert get_marker(StringIO(no_hits), orf_finder='prodigal', multi_mode=True) == {}
    rs = get_marker(StringIO(no_hits), orf_finder='prodigal', contig_to_marker=True)
    assert rs['k119_224042'] == []


def test_get_marker_multiple():
    rs = get_marker('test/data/concatenated.hmmout.gz', multi_mode=True, orf_finder='fraggenescan')
    rs
    # This was computed by running the previous version
    assert rs == {
         'bin000005': ['k119_46822'],
         'bin000051': ['k119_368830', 'k119_419149'],
         'bin000052': ['k119_314249', 'k119_686387'],
         'bin000056': ['k119_266899', 'k119_447637'],
         'bin000057': ['k119_464254'],
         'bin000078': ['k119_162461', 'k119_278484', 'k119_678785', 'k119_684640'],
         'bin000085': ['k119_232333', 'k119_549069'],
         'bin000110': ['k119_276145', 'k119_50857']}

def test_estimate_seeds():
    seeds = estimate_seeds('test/single_sample_data/input.fasta',
                                1000, num_process=1, orf_finder='fast-naive')
    assert len(seeds) > 0


def _run_estimate_seeds(out, **kwargs):
    defaults = dict(
        fasta_path='test/single_sample_data/input.fasta',
        binned_length=1000,
        num_process=1,
        orf_finder='fast-naive',
        output=out,
    )
    defaults.update(kwargs)
    return estimate_seeds(**defaults)


def test_marker_cache_hit_skips_hmmsearch(tmp_path):
    out = str(tmp_path / 'out')
    _run_estimate_seeds(out)
    assert os.path.exists(os.path.join(out, 'markers.hmmout'))
    assert os.path.exists(os.path.join(out, 'markers.hmmout.json'))

    with mock.patch('SemiBin.markers.subprocess.check_call') as check_call:
        _run_estimate_seeds(out)
        assert check_call.call_count == 0


def test_marker_cache_invalidates_on_param_change(tmp_path):
    out = str(tmp_path / 'out')
    _run_estimate_seeds(out, binned_length=1000)
    sidecar = os.path.join(out, 'markers.hmmout.json')
    with open(sidecar) as f:
        first = json.load(f)
    assert first['binned_length'] == 1000

    _run_estimate_seeds(out, binned_length=1500)
    with open(sidecar) as f:
        second = json.load(f)
    assert second['binned_length'] == 1500


def test_marker_cache_invalidates_when_sidecar_missing(tmp_path):
    out = str(tmp_path / 'out')
    _run_estimate_seeds(out)
    os.remove(os.path.join(out, 'markers.hmmout.json'))

    with mock.patch('SemiBin.markers.subprocess.check_call',
                    wraps=__import__('subprocess').check_call) as check_call:
        _run_estimate_seeds(out)
        assert check_call.call_count == 1
    assert os.path.exists(os.path.join(out, 'markers.hmmout.json'))


def test_marker_cache_invalidates_on_fasta_change(tmp_path):
    fasta = tmp_path / 'input.fasta'
    shutil.copy('test/single_sample_data/input.fasta', fasta)
    out = str(tmp_path / 'out')
    _run_estimate_seeds(out, fasta_path=str(fasta))

    # Bump mtime to a clearly different value to simulate a regenerated FASTA.
    st = os.stat(fasta)
    os.utime(fasta, ns=(st.st_atime_ns, st.st_mtime_ns + 10_000_000_000))

    with mock.patch('SemiBin.markers.subprocess.check_call',
                    wraps=__import__('subprocess').check_call) as check_call:
        _run_estimate_seeds(out, fasta_path=str(fasta))
        assert check_call.call_count == 1


def test_n_single_copy_markers_matches_hmm():
    with open(markers.MARKER_HMM_PATH) as f:
        names = [line.split()[1] for line in f if line.startswith('NAME ')]
    normalize = markers.__normalize_marker_trans
    assert set(normalize.keys()) <= set(names)
    assert set(normalize.values()) <= set(names)
    assert len({normalize.get(n, n) for n in names}) == markers.N_SINGLE_COPY_MARKERS
