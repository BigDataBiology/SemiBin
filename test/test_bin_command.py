import os
import pandas as pd


### Input fa
os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --model test/bin_data/model.h5 -i test/bin_data/input.fasta -o output_bin_fa -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_fa/output_bins')) > 0
assert len(os.listdir('output_bin_fa/output_recluster_bins')) > 0


### Input .gz
os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --model test/bin_data/model.h5 -i test/bin_data/input.fasta.gz -o output_bin_gz -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_gz/output_bins')) > 0
assert len(os.listdir('output_bin_gz/output_recluster_bins')) > 0

### Input .bz2
os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --model test/bin_data/model.h5 -i test/bin_data/input.fasta.bz2 -o output_bin_bz2 -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_bz2/output_bins')) > 0
assert len(os.listdir('output_bin_bz2/output_recluster_bins')) > 0

### Input .xz
os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --model test/bin_data/model.h5 -i test/bin_data/input.fasta.xz -o output_bin_xz -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_xz/output_bins')) > 0
assert len(os.listdir('output_bin_xz/output_recluster_bins')) > 0

### use pretrained model
os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --environment human_gut -i test/bin_data/input.fasta.xz -o output_bin_human_gut -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_human_gut/output_bins')) > 0
assert len(os.listdir('output_bin_human_gut/output_recluster_bins')) > 0

os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --environment ocean -i test/bin_data/input.fasta.xz -o output_bin_ocean -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_ocean/output_bins')) > 0
assert len(os.listdir('output_bin_ocean/output_recluster_bins')) > 0

os.system('SemiBin bin --data test/bin_data/data.csv --minfasta-kbs 200 --recluster --max-edges 20 --max-node 1 --environment dog_gut -i test/bin_data/input.fasta.xz -o output_bin_dog_gut -m 2500 --ratio 0.05 -p 1')

assert len(os.listdir('output_bin_dog_gut/output_bins')) > 0
assert len(os.listdir('output_bin_dog_gut/output_recluster_bins')) > 0
