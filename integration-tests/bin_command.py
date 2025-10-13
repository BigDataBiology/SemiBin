import os
from os import path
import subprocess
import pandas as pd
import shutil
from glob import glob
from SemiBin.fasta import fasta_iter


def check_bin_output(fafile, odir):
    contigs_original = dict(fasta_iter(fafile))
    bin_info = pd.read_table(f'{odir}/contig_bins.tsv')
    contig2bin = bin_info.set_index('contig').to_dict()['bin']
    if path.exists(f'{odir}/bins_info.tsv'):
        tab = pd.read_table(f'{odir}/bins_info.tsv', index_col=0)
    else:
        tab = pd.read_table(f'{odir}/recluster_bins_info.tsv', index_col=0)
    for f in tab.index:
        assert 'output_bins' in f
        assert os.path.exists(f)

    contigs_binned = []
    sequences_binned = {}
    n_bins = 0
    for bf in os.listdir(f'{odir}/output_bins/'):
        if not bf.endswith('.gz'):
            continue
        n_bins += 1
        bin_id = None
        for h, seq in fasta_iter(f'{odir}/output_bins/{bf}'):
            if bin_id is None:
                bin_id = contig2bin[h]
            else:
                assert bin_id == contig2bin[h]
            contigs_binned.append(h)
            sequences_binned[h] = seq
    assert len(set(bin_info['bin']) - set([-1])) == n_bins

    # Check that contigs are not repeated across bins (or within a bin)
    assert len(contigs_binned) == len(set(contigs_binned))

    # Check that contigs do not change (and no new contigs are added)
    for contig in sequences_binned:
        assert contigs_original[contig] == sequences_binned[contig]


ifile = 'input.fasta.xz'
odir = 'test-outputs/output_bin_xz'
shutil.rmtree(odir, ignore_errors=True)
subprocess.check_call(
    ['SemiBin2', 'bin_short',
     '--data', 'test/bin_data/data.csv',
     '--minfasta-kbs', '0',
     '--max-edges', '20',
     '--max-node', '1',
     '--model', 'test/bin_data/model.pt',
     '-i', f'test/bin_data/{ifile}',
     '-o', odir,
     '-m', '2500',
     '--ratio', '0.05',])

assert len(os.listdir(f'{odir}/output_bins')) > 0
check_bin_output(f'test/bin_data/{ifile}', odir)


odir = f'test-outputs/output_long'
shutil.rmtree(odir, ignore_errors=True)
subprocess.check_call(
    ['SemiBin2', 'bin_long',
     '--data', 'test/bin_data/data.csv',
     '--minfasta-kbs', '0',
     '--model', 'test/bin_data/model.pt',
     '-i', 'test/bin_data/input.fasta',
     '-o', odir,
     '-m', '2500',
     '--ratio', '0.05',
     '--tag-output', 'TestTag'
     ])
assert len(glob(f'{odir}/output_bins/TestTag_*')) > 0
check_bin_output(f'test/bin_data/input.fasta', odir)


# Different pretrained models
for env,odir in [
        ('human_gut', 'output_human_gut'),
        ('dog_gut', 'output_dog_gut'),
        ('ocean', 'output_ocean'),
        ]:
    odir = f'test-outputs/{odir}'
    shutil.rmtree(odir, ignore_errors=True)
    subprocess.check_call(
        ['SemiBin2', 'bin_short',
         '--data', 'test/bin_data/data.csv',
         '--minfasta-kbs', '0',
         '--max-edges', '20',
         '--max-node', '1',
         '--environment', env,
         '-i', 'test/bin_data/input.fasta',
         '-o', odir,
         '-m', '2500',
         '--ratio', '0.05'])
    assert len(os.listdir(odir+'/output_bins')) > 0
    check_bin_output('test/bin_data/input.fasta', odir)


