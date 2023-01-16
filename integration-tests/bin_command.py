import os
import subprocess
import pandas as pd
import shutil
from glob import glob

# Test different input formats
for ifile, odir in [
        ('input.fasta', 'output_bin_fa'),
        ('input.fasta.gz', 'output_bin_gz'),
        ('input.fasta.bz2', 'output_bin_bz2'),
        ('input.fasta.xz', 'output_bin_xz'),
        ]:
    odir = f'test-outputs/{odir}_no_recluster'
    shutil.rmtree(odir, ignore_errors=True)
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
         '--no-recluster',
         '-p', '1'])

    assert len(os.listdir(f'{odir}/output_bins')) > 0

odir = f'test-outputs/output_long'
shutil.rmtree(odir, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'bin_long',
     '--data', 'test/bin_data/data.csv',
     '--minfasta-kbs', '0',
     '--model', 'test/bin_data/model.h5',
     '-i', 'test/bin_data/input.fasta',
     '-o', odir,
     '-m', '2500',
     '--ratio', '0.05',
     '--tag-output', 'TestTag',
     '-p', '1'])
assert len(glob(f'{odir}/output_bins/TestTag_*')) > 0


ifile = 'input.fasta'
odir = 'test-outputs/with_recluster'
shutil.rmtree(odir, ignore_errors=True)
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

assert len(os.listdir(f'{odir}/output_prerecluster_bins')) > 0
assert len(os.listdir(f'{odir}/output_recluster_bins')) > 0

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
    assert len(os.listdir(odir+'/output_prerecluster_bins')) > 0
    assert len(os.listdir(odir+'/output_recluster_bins')) > 0
    tab = pd.read_table(f'{odir}/recluster_bins_info.tsv', index_col=0)
    for f in tab.index:
        assert 'output_recluster_bins' in f
        assert os.path.exists(f)


# Test with input taxonomy file

single_sample_input = 'test/single_sample_data'
multi_sample_input =  'test/multi_samples_data'
single_cannot_output = 'test-outputs/single_output_cannot_with_taxonomy'
single_cannot_ref_output = 'test-outputs/single_output_cannot_with_ref'

shutil.rmtree(single_cannot_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'generate_cannot_links',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_cannot_output,
     '--taxonomy-annotation-table', f'{single_sample_input}/taxonomyResult.tsv'])
assert os.path.exists(f'{single_cannot_output}/cannot/cannot.txt')

shutil.rmtree(single_cannot_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'generate_cannot_links',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_cannot_ref_output,
     f'-r{single_sample_input}/reference_genome'])
assert os.path.exists(f'{single_cannot_ref_output}/cannot/cannot.txt')
