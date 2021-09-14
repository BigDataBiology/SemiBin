import os
import pandas as pd


### Input fa
os.system('SemiBin generate_data_single -i test/single_sample_data/input.fasta -o output_single_fa -m 2500 --ratio 0.05 --ml-threshold 4000 -p 1 -b test/single_sample_data/input.sorted.bam')

data = pd.read_csv('output_single_fa/data.csv', index_col=0)
data_split = pd.read_csv('output_single_fa/data_split.csv', index_col=0)

assert data.shape == (40, 138)
assert data_split.shape == (80, 136)
