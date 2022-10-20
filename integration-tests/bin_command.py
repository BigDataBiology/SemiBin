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
    assert len(os.listdir(f'{odir}/output_prerecluster_bins')) > 0
    assert len(os.listdir(f'{odir}/output_recluster_bins')) > 0

    odir = f'test-outputs/{odir}_long'
    subprocess.check_call(
        ['SemiBin', 'bin_long',
         '--data', 'test/bin_data/data.csv',
         '--minfasta-kbs', '0',
         '--model', 'test/bin_data/model.h5',
         '-i', f'test/bin_data/{ifile}',
         '-o', odir,
         '-m', '2500',
         '--ratio', '0.05',
         '-p', '1'])
    assert len(os.listdir(f'{odir}/output_bins')) > 0


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
    assert len(os.listdir(odir+'/output_prerecluster_bins')) > 0
    assert len(os.listdir(odir+'/output_recluster_bins')) > 0

# Test with input taxonomy file

single_sample_input = 'test/single_sample_data'
multi_sample_input =  'test/multi_samples_data'
single_cannot_output = 'test-outputs/single_output_cannot_with_taxonomy'
single_cannot_ref_output = 'test-outputs/single_output_cannot_with_ref'
single_output = 'test-outputs/single_output_with_taxonomy'
single_output_ref = 'test-outputs/single_output_with_ref'
multi_output = 'test-outputs/multi_output_with_taxonomy'
multi_output_ref = 'test-outputs/multi_output_with_ref'

subprocess.check_call(f'SemiBin generate_cannot_links -i {single_sample_input}/input.fasta -o {single_cannot_output} --taxonomy-annotation-table {single_sample_input}/taxonomyResult.tsv', shell=True)
assert os.path.exists(f'{single_cannot_output}/cannot/cannot.txt')

subprocess.check_call(f'SemiBin generate_cannot_links -i {single_sample_input}/input.fasta -o {single_cannot_ref_output} -r{single_sample_input}/reference_genome ', shell=True)
assert os.path.exists(f'{single_cannot_ref_output}/cannot/cannot.txt')

subprocess.check_call(f'SemiBin single_easy_bin -i {single_sample_input}/input.fasta -o {single_output} -b {single_sample_input}/input.sorted.bam --taxonomy-annotation-table {single_sample_input}/taxonomyResult.tsv --epoches 1 --training-type semi --sequencing-type short_read', shell=True)
assert os.path.exists(f'{single_output}/output_prerecluster_bins')
assert os.path.exists(f'{single_output}/output_recluster_bins')

subprocess.check_call(f'SemiBin single_easy_bin -i {single_sample_input}/input.fasta -o {single_output_ref} -b {single_sample_input}/input.sorted.bam -r {single_sample_input}/reference_genome --epoches 1 --training-type semi --sequencing-type short_read', shell=True)
assert os.path.exists(f'{single_output_ref}/output_prerecluster_bins')
assert os.path.exists(f'{single_output_ref}/output_recluster_bins')

subprocess.check_call(f'SemiBin multi_easy_bin -i {multi_sample_input}/input_multi.fasta -o {multi_output} -b {multi_sample_input}/*.bam --taxonomy-annotation-table {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv {single_sample_input}/taxonomyResult.tsv -s : --epoches 1 --training-type semi --sequencing-type short_read',  shell=True)
assert os.path.exists(f'{multi_output}/bins')
for i in range(10):
    assert os.path.exists(f'{multi_output}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_output}/samples/S{i + 1}/output_recluster_bins')


subprocess.check_call(f'SemiBin multi_easy_bin -i {multi_sample_input}/input_multi.fasta -o {multi_output_ref} -b {multi_sample_input}/*.bam -r {single_sample_input}/reference_genome -s : --epoches 1 --training-type semi --sequencing-type short_read',  shell=True)
assert os.path.exists(f'{multi_output_ref}/bins')
for i in range(10):
    assert os.path.exists(f'{multi_output_ref}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_output_ref}/samples/S{i+1}/output_recluster_bins')


# Test .cram format input
single_cram_output = 'test-outputs/single_output_cram'
subprocess.check_call(f'SemiBin single_easy_bin -i {single_sample_input}/input.fasta -o {single_cram_output} -b {single_sample_input}/input.cram --environment human_gut --training-type semi --sequencing-type short_read', shell=True)
assert os.path.exists(f'{single_cram_output}/output_prerecluster_bins')
assert os.path.exists(f'{single_cram_output}/output_recluster_bins')

multi_output_cram = 'test-outputs/multi_output_cram'
subprocess.check_call(f'SemiBin multi_easy_bin -i {multi_sample_input}/input_multi.fasta -o {multi_output_cram} -b {multi_sample_input}/*.cram -r {single_sample_input}/reference_genome -s : --epoches 1 --training-type semi --sequencing-type short_read',  shell=True)
assert os.path.exists(f'{multi_output_cram}/bins')
for i in range(10):
    assert os.path.exists(f'{multi_output_cram}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_output_cram}/samples/S{i+1}/output_recluster_bins')

# Test training with self-supervised learning
single_self_output = 'test-outputs/single_output_self'
subprocess.check_call(f'SemiBin single_easy_bin -i {single_sample_input}/input.fasta -o {single_self_output} -b {single_sample_input}/input.sorted.bam --epoches 1 --training-type self --sequencing-type short_read', shell=True)
assert os.path.exists(f'{single_self_output}/output_prerecluster_bins')
assert os.path.exists(f'{single_self_output}/output_recluster_bins')

multi_self_output = 'test-outputs/multi_output_self'
subprocess.check_call(f'SemiBin multi_easy_bin -i {multi_sample_input}/input_multi.fasta -o {multi_self_output} -b {multi_sample_input}/*.bam -s : --epoches 1 --training-type self --sequencing-type short_read',  shell=True)
assert os.path.exists(f'{multi_self_output}/bins')
for i in range(10):
    assert os.path.exists(f'{multi_self_output}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_self_output}/samples/S{i+1}/output_recluster_bins')

# Test binning with long-read
single_self_output_long = 'test-outputs/single_output_self_long'
subprocess.check_call(f'SemiBin single_easy_bin -i {single_sample_input}/input.fasta -o {single_self_output_long} -b {single_sample_input}/input.sorted.bam --epoches 1 --training-type self --sequencing-type long_read', shell=True)
assert os.path.exists(f'{single_self_output_long}/output_bins')

multi_self_output_long = 'test-outputs/multi_output_self_long'
subprocess.check_call(f'SemiBin multi_easy_bin -i {multi_sample_input}/input_multi.fasta -o {multi_self_output_long} -b {multi_sample_input}/*.bam -s : --epoches 1 --training-type self --sequencing-type long_read',  shell=True)
assert os.path.exists(f'{multi_self_output_long}/bins')
for i in range(10):
    assert os.path.exists(f'{multi_self_output_long}/samples/S{i+1}/output_bins')