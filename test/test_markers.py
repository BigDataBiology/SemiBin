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
