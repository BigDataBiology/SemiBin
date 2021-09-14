import os
import pandas as pd


### Input fa
os.system('SemiBin generate_data_single -i test/coassembly_sample_data/input.fasta -o output_coassembly_fa -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/coassembly_sample_data/input.sorted*.bam')

data = pd.read_csv('output_coassembly_fa/data.csv', index_col=0)
data_split = pd.read_csv('output_coassembly_fa/data_split.csv', index_col=0)

assert data.shape == (40, 138)
assert data_split.shape == (80, 136)


### Input .gz
os.system('SemiBin generate_data_single -i test/coassembly_sample_data/input.fasta.gz -o output_coassembly_gz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/coassembly_sample_data/input.sorted*.bam')

data = pd.read_csv('output_coassembly_gz/data.csv', index_col=0)
data_split = pd.read_csv('output_coassembly_gz/data_split.csv', index_col=0)

assert data.shape == (40, 138)
assert data_split.shape == (80, 136)

### Input .bz2
os.system('SemiBin generate_data_single -i test/coassembly_sample_data/input.fasta.bz2 -o output_coassembly_bz2 -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/coassembly_sample_data/input.sorted*.bam')

data = pd.read_csv('output_coassembly_bz2/data.csv', index_col=0)
data_split = pd.read_csv('output_coassembly_bz2/data_split.csv', index_col=0)

assert data.shape == (40, 138)
assert data_split.shape == (80, 136)

### Input .xz
os.system('SemiBin generate_data_single -i test/coassembly_sample_data/input.fasta.xz -o output_coassembly_xz -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/coassembly_sample_data/input.sorted*.bam')

data = pd.read_csv('output_coassembly_xz/data.csv', index_col=0)
data_split = pd.read_csv('output_coassembly_xz/data_split.csv', index_col=0)

assert data.shape == (40, 138)
assert data_split.shape == (80, 136)