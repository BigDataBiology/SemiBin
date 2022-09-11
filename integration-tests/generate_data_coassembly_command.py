import os
import pandas as pd
import subprocess
from glob import glob

os.makedirs('test-outputs', exist_ok=True)
bam_files = glob('test/coassembly_sample_data/input.sorted*.bam')


for ext in ['', '.gz', '.bz2', '.xz']:
    odir = f"test-outputs/output_coassembly{'_'+ext.replace('.', '_') if ext else '_fa'}"
    subprocess.check_call(
        ['SemiBin', 'generate_sequence_features_single',
         '-i', f'test/coassembly_sample_data/input.fasta{ext}',
         '-o', odir,
         '-m', '2500',
         '--ratio', '0.05',
         '--ml-threshold', '4000',
         '-p', '1',
         '-b'] + bam_files)

    data = pd.read_csv(f'{odir}/data.csv', index_col=0)
    data_split = pd.read_csv(f'{odir}/data_split.csv', index_col=0)

    assert data.shape == (40, 141)
    assert data_split.shape == (80, 141)


