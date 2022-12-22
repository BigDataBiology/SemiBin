import os
import subprocess
import pandas as pd
import shutil

def sglob(pat):
    from glob import glob
    r = glob(pat)
    r.sort()
    return r

single_sample_input = 'test/single_sample_data'
multi_sample_input =  'test/multi_samples_data'

single_output = 'test-outputs/single_output_with_taxonomy'
single_output_ref = 'test-outputs/single_output_with_ref'
multi_output = 'test-outputs/multi_output_with_taxonomy'
multi_output_ref = 'test-outputs/multi_output_with_ref'

shutil.rmtree(single_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'single_easy_bin',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_output,
     '-b', f'{single_sample_input}/input.sorted.bam',
     '--taxonomy-annotation-table', f'{single_sample_input}/taxonomyResult.tsv',
     '--epochs', '1',
     '--semi-supervised'])
assert os.path.exists(f'{single_output}/output_prerecluster_bins')
assert os.path.exists(f'{single_output}/output_recluster_bins')

shutil.rmtree(single_output_ref, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'single_easy_bin',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_output_ref,
     '-b', f'{single_sample_input}/input.sorted.bam',
     '-r', f'{single_sample_input}/reference_genome',
     '--epochs', '1'])
assert os.path.exists(f'{single_output_ref}/output_prerecluster_bins')
assert os.path.exists(f'{single_output_ref}/output_recluster_bins')

multi_sample_input_bams = sglob(f'{multi_sample_input}/*.bam')
shutil.rmtree(multi_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin',
     'multi_easy_bin',
     '-i', f'{multi_sample_input}/input_multi.fasta',
     '-o', multi_output,
     '-b'] + multi_sample_input_bams + [
     '--taxonomy-annotation-table'] + [
         f'{single_sample_input}/taxonomyResult.tsv' for _ in multi_sample_input_bams] + [
     '-s', ':',
     '--epochs', '1',
     '--semi-supervised'])

assert os.path.exists(f'{multi_output}/bins')
for i in range(len(multi_sample_input_bams)):
    assert os.path.exists(f'{multi_output}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_output}/samples/S{i+1}/output_recluster_bins')

shutil.rmtree(multi_output_ref, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'multi_easy_bin',
     '-i', f'{multi_sample_input}/input_multi.fasta',
     '-o', multi_output_ref,
     '-b'] + multi_sample_input_bams + [
     '-r', f'{single_sample_input}/reference_genome',
     '-s', ':',
     '--epochs', '1',
     '--semi-supervised'])
assert os.path.exists(f'{multi_output_ref}/bins')
for i in range(len(multi_sample_input_bams)):
    assert os.path.exists(f'{multi_output_ref}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_output_ref}/samples/S{i+1}/output_recluster_bins')

# Test .cram format input
single_cram_output = 'test-outputs/single_output_cram'
shutil.rmtree(single_cram_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'single_easy_bin',
     '-i', f'{single_sample_input}/input.fasta',
     '-b', f'{single_sample_input}/input.cram',
     '-o', single_cram_output,
     '--environment', 'human_gut'
     ])
assert os.path.exists(f'{single_cram_output}/output_prerecluster_bins')

# Test binning with long-read
single_self_output_long = 'test-outputs/single_output_self_long'
shutil.rmtree(single_self_output_long, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'single_easy_bin',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_self_output_long,
     '-b', f'{single_sample_input}/input.sorted.bam',
     '--epochs', '1',
     '--self-supervised',
     '--orf-finder', 'fast-naive',
     '--sequencing-type', 'long_read'])
assert os.path.exists(f'{single_self_output_long}/output_bins')

multi_self_output_long = 'test-outputs/multi_output_self_long'
shutil.rmtree(multi_self_output_long, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'multi_easy_bin',
     '-i', f'{multi_sample_input}/input_multi.fasta',
     '-o', multi_self_output_long,
     '-b'] + multi_sample_input_bams + [
     '-s', ':',
     '--epochs', '1',
     '--training-type', 'self',
     '--sequencing-type', 'long_read'])
assert os.path.exists(f'{multi_self_output_long}/bins')
for i in range(len(multi_sample_input_bams)):
    assert os.path.exists(f'{multi_self_output_long}/samples/S{i+1}/output_bins')
assert os.path.exists(f'{single_cram_output}/output_recluster_bins')


# Test training with self-supervised learning
single_self_output = 'test-outputs/single_output_self'
shutil.rmtree(single_self_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin', 'single_easy_bin',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_self_output,
     '-b', f'{single_sample_input}/input.sorted.bam',
     '--epochs', '1',
     '--self-supervised'])
assert os.path.exists(f'{single_self_output}/output_prerecluster_bins')
assert os.path.exists(f'{single_self_output}/output_recluster_bins')

multi_self_output = 'test-outputs/multi_output_self'
shutil.rmtree(multi_self_output, ignore_errors=True)
subprocess.check_call(
        ['SemiBin', 'multi_easy_bin',
        '-i', f'{multi_sample_input}/input_multi.fasta',
        '-o', multi_self_output,
        '-b'] + multi_sample_input_bams + [
        '-s', ':',
        '--epochs', '1',
        '--self-supervised'])

assert os.path.exists(f'{multi_self_output}/bins')
for i in range(len(multi_sample_input_bams)):
    assert os.path.exists(f'{multi_self_output}/samples/S{i+1}/output_prerecluster_bins')
    assert os.path.exists(f'{multi_self_output}/samples/S{i+1}/output_recluster_bins')
