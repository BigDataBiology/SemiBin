import os
import subprocess
import pandas as pd
import shutil
from glob import glob
from SemiBin.fasta import fasta_iter

def test_overlap(contig, bin_dir, prefix):
    contig_dict = {}
    for h, seq in fasta_iter(contig):
        contig_dict[h] = seq

    contig_bin_name = []
    contig_bin_dict = {}
    bin_list = os.listdir(bin_dir)
    bin_list = [temp for temp in bin_list if temp.split('.')[-1] == '{}'.format(prefix)]

    for bin in bin_list:
        for h, seq in fasta_iter(
                '{0}/{1}'.format(bin_dir, bin)):
            contig_bin_name.append(h)
            contig_bin_dict[h] = seq

    assert len(contig_bin_name) == len(set(contig_bin_name))

    for contig in contig_bin_dict:
        assert contig_dict[contig] == contig_bin_dict[contig]

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
        ['SemiBin2', 'bin_short',
         '--data', 'test/bin_data/data.csv',
         '--minfasta-kbs', '0',
         '--max-edges', '20',
         '--max-node', '1',
         '--model', 'test/bin_data/model.pt',
         '-i', f'test/bin_data/{ifile}',
         '-o', odir,
         '-m', '2500',
         '--ratio', '0.05',
         '-p', '1'])

    assert len(os.listdir(f'{odir}/output_bins')) > 0
    test_overlap(f'test/bin_data/{ifile}', f'{odir}/output_bins', "gz")

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
     '--tag-output', 'TestTag',
     '-p', '1'])
assert len(glob(f'{odir}/output_bins/TestTag_*')) > 0
test_overlap(f'test/bin_data/input.fasta', f'{odir}/output_bins', "gz")


# Different pretrained models
for env,odir in [
        ('human_gut', 'output_human_gut'),
        ('dog_gut', 'output_dog_gut'),
        ('ocean', 'output_ocean'),
        ]:
    odir = f'test-outputs/{odir}'
    subprocess.check_call(
        ['SemiBin2', 'bin_short',
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
    tab = pd.read_table(f'{odir}/recluster_bins_info.tsv', index_col=0)
    for f in tab.index:
        assert 'output_bins' in f
        assert os.path.exists(f)


