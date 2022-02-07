import subprocess
import pandas as pd

### Input fa
subprocess.check_call('SemiBin generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta -o output_multi_fa -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('output_multi_fa/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'output_multi_fa/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)


### Input .gz
subprocess.check_call('SemiBin generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.gz -o output_multi_gz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('output_multi_gz/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'output_multi_gz/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)

### Input .bz2
subprocess.check_call('SemiBin generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.bz2 -o output_multi_bz2 -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('output_multi_bz2/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'output_multi_bz2/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)

### Input .xz
subprocess.check_call('SemiBin generate_sequence_features_multi -i test/multi_samples_data/input_multi.fasta.xz -o output_multi_xz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/multi_samples_data/input_multi_sorted*.bam -s :', shell=True)

for i in range(10):
    data = pd.read_csv('output_multi_xz/samples/S{}/data.csv'.format(i + 1),
                       index_col=0)
    data_split = pd.read_csv(
        'output_multi_xz/samples/S{}/data_split.csv'.format(i + 1), index_col=0)
    assert data.shape == (20, 146)
    assert data_split.shape == (40, 146)
