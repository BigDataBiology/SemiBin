import os
import subprocess

# Test different input formats
for ifile, odir in [
        ('input.fasta', 'output_bin_fa'),
        ('input.fasta.gz', 'output_bin_gz'),
        ('input.fasta.bz2', 'output_bin_bz2'),
        ('input.fasta.xz', 'output_bin_xz'),
        ]:
    odir = f'test-outputs/{odir}'
    subprocess.check_call(
        ['SemiBin', 'bin',
         '--data', 'test/bin_data/data.csv',
         '--minfasta-kbs', '0',
         '--max-edges', '20',
         '--max-node', '1',
         '--model', 'test/bin_data/model.h5',
         '-i', f'test/bin_data/{ifile}',
         '-o', odir,
         '-m', '2500',
         '--ratio', '0.05',
         '-p', '1'])
    assert len(os.listdir(f'{odir}/output_bins')) > 0
    assert len(os.listdir(f'{odir}/output_recluster_bins')) > 0

ifile = 'input.fasta'
odir = 'test-outputs/no_recluster'
subprocess.check_call(
    ['SemiBin', 'bin',
     '--data', 'test/bin_data/data.csv',
     '--minfasta-kbs', '0',
     '--max-edges', '20',
     '--max-node', '1',
     '--no-recluster',
     '--model', 'test/bin_data/model.h5',
     '-i', f'test/bin_data/{ifile}',
     '-o', odir,
     '-m', '2500',
     '--ratio', '0.05',
     '-p', '1'])
assert len(os.listdir(f'{odir}/output_bins')) > 0
assert not os.path.exists(f'{odir}/output_recluster_bins')

# Different pretrained models
for env,odir in [
        ('human_gut', 'output_human_gut'),
        ('dog_gut', 'output_dog_gut'),
        ('ocean', 'output_ocean'),
        ]:
    odir = f'test-outputs/{odir}'
    subprocess.check_call(
        ['SemiBin', 'bin',
         '--data', 'test/bin_data/data.csv',
         '--minfasta-kbs', '0',
         '--max-edges', '20',
         '--max-node', '1',
         '--environment', env,
         '-i', 'test/bin_data/input.fasta.xz',
         '-o', odir,
         '-m', '2500',
         '--ratio', '0.05',
         '-p', '1'])
    assert len(os.listdir(odir+'/output_bins')) > 0
    assert len(os.listdir(odir+'/output_recluster_bins')) > 0
