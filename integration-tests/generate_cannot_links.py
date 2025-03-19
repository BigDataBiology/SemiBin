import subprocess
import shutil
import os

# Test with input taxonomy file

single_sample_input = 'test/single_sample_data'
multi_sample_input =  'test/multi_samples_data'
single_cannot_output = 'test-outputs/single_output_cannot_with_taxonomy'
single_cannot_ref_output = 'test-outputs/single_output_cannot_with_ref'

shutil.rmtree(single_cannot_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin2', 'generate_cannot_links',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_cannot_output,
     '--taxonomy-annotation-table', f'{single_sample_input}/taxonomyResult.tsv'])
assert os.path.exists(f'{single_cannot_output}/cannot/cannot.txt')

shutil.rmtree(single_cannot_output, ignore_errors=True)
subprocess.check_call(
    ['SemiBin2', 'generate_cannot_links',
     '-i', f'{single_sample_input}/input.fasta',
     '-o', single_cannot_ref_output,
     f'-r{single_sample_input}/reference_genome'])
assert os.path.exists(f'{single_cannot_ref_output}/cannot/cannot.txt')
