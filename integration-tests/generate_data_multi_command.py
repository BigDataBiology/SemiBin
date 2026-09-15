import subprocess
import pandas as pd

# NOTE: test/multi_samples_data/input_multi.fasta (and its .gz/.bz2/.xz variants)
# only contains 3 samples (S1-S3) to keep these tests fast. The BAM/abundance
# fixtures still contain 10 samples' worth of coverage data (used to keep the
# is_combined/n_sample>=5 code paths exercised); generate_sequence_features_multi
# only outputs data for the sample groups present in the input FASTA, so the
# extra coverage data for S4-S10 is simply ignored.
# test/multi_samples_data/input_multi_full.fasta is the original untruncated
# (10-sample) fixture, kept around only for the CRAM test below, since the
# .cram files were encoded against the full reference and need all of its
# contigs to decode.
N_SAMPLES = 3

### Input fa
subprocess.check_call('SemiBin2 generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta -o test-outputs/output_multi_fa -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(N_SAMPLES):
    data = pd.read_csv('test-outputs/output_multi_fa/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_fa/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)

# running with abundance file from strobealign-aemb
subprocess.check_call('SemiBin2 generate_sequence_features_multi '
                      '-i test/multi_samples_data/input_multi.fasta '
                      '-o test-outputs/output_multi_fa -m 2500 '
                      '--ratio 0.05 --ml-threshold 4000 -p 1 '
                      '-a test/multi_samples_data/*.txt -s :', shell=True)

for i in range(N_SAMPLES):
    data = pd.read_csv('test-outputs/output_multi_fa/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_fa/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)



### Input .gz (compressed-input handling is otherwise identical to .fa above, so just
### smoke-test that reading a compressed multi-sample FASTA works; .bz2/.xz are covered
### by generate_data_single_command.py, which already exercises those code paths.)
subprocess.check_call('SemiBin2 generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.gz -o test-outputs/output_multi_gz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(N_SAMPLES):
    data = pd.read_csv('test-outputs/output_multi_gz/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_gz/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)


multi_sample_input =  'test/multi_samples_data'
multi_output_cram = 'test-outputs/multi_output_cram'
# Uses the untruncated (10-sample) fixture: the .cram files were encoded against
# it, so all of its contigs are needed to decode them.
subprocess.check_call(f'SemiBin2 generate_sequence_features_multi -i {multi_sample_input}/input_multi_full.fasta -o {multi_output_cram} -b {multi_sample_input}/*.cram -s :', shell=True)

for i in range(10):
    data = pd.read_csv(
        f'{multi_output_cram}/samples/S{i+1}/data.csv', index_col=0)
    data_split = pd.read_csv(
        f'{multi_output_cram}/samples/S{i+1}/data_split.csv', index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)
