import subprocess
import pandas as pd

### Input fa
subprocess.check_call('SemiBin2 generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta -o test-outputs/output_multi_fa -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
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

for i in range(10):
    data = pd.read_csv('test-outputs/output_multi_fa/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_fa/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)



### Input .gz
subprocess.check_call('SemiBin2 generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.gz -o test-outputs/output_multi_gz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('test-outputs/output_multi_gz/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_gz/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)

### Input .bz2
subprocess.check_call('SemiBin2 generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.bz2 -o test-outputs/output_multi_bz2 -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('test-outputs/output_multi_bz2/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_bz2/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)

### Input .xz
subprocess.check_call('SemiBin2 generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.xz -o test-outputs/output_multi_xz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('test-outputs/output_multi_xz/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'test-outputs/output_multi_xz/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)


multi_sample_input =  'test/multi_samples_data'
multi_output_cram = 'test-outputs/multi_output_cram'
subprocess.check_call(f'SemiBin2 generate_sequence_features_multi -i {multi_sample_input}/input_multi.fasta -o {multi_output_cram} -b {multi_sample_input}/*.cram -s :', shell=True)

for i in range(10):
    data = pd.read_csv(
        f'{multi_output_cram}/samples/S{i+1}/data.csv', index_col=0)
    data_split = pd.read_csv(
        f'{multi_output_cram}/samples/S{i+1}/data_split.csv', index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)

