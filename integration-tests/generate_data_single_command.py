import pandas as pd
import subprocess

data0 = None
data0_split = None
for ifile in [
        'input.fasta',
        'input.fasta.gz',
        'input.fasta.bz2',
        'input.fasta.xz',
        ]:
    odir_suf = ifile.split('.')[1]
    subprocess.check_call(
        ['SemiBin2', 'generate_sequence_features_single',
         '-i', f'test/single_sample_data/{ifile}',
         '-o', f'test-outputs/output_single_{odir_suf}',
         '-m', '2500',
         '--ratio', '0.05',
         '--ml-threshold', '4000',
         '-p', '1',
         '-b', 'test/single_sample_data/input.sorted.bam'])

    data = pd.read_csv(f'test-outputs/output_single_{odir_suf}/data.csv', index_col=0)
    data_split = pd.read_csv(f'test-outputs/output_single_{odir_suf}/data_split.csv', index_col=0)

    assert data.shape == (40, 138)
    assert data_split.shape == (80, 136)
    if data0 is not None:
        assert (data == data0).all().all()
        assert (data_split == data0_split).all().all()
    else:
        data0 = data
        data0_split = data_split

