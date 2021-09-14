import os
import pandas as pd


### Input fa
os.system('SemiBin train --data test/train_data/data.csv --data-split test/train_data/data_split.csv -c test/train_data/cannot.txt --epoches 1 --batch-size 2048 --mode single -i test/train_data/input.fasta -o output_train_fa -m 2500 --ratio 0.05 -p 1')

assert os.path.exists('output_train_fa/model.h5')

### Input .gz
os.system('SemiBin train --data test/train_data/data.csv --data-split test/train_data/data_split.csv -c test/train_data/cannot.txt --epoches 1 --batch-size 2048 --mode single -i test/train_data/input.fasta.gz -o output_train_gz -m 2500 --ratio 0.05 -p 1')

assert os.path.exists('output_train_gz/model.h5')

### Input .bz2
os.system('SemiBin train --data test/train_data/data.csv --data-split test/train_data/data_split.csv -c test/train_data/cannot.txt --epoches 1 --batch-size 2048 --mode single -i test/train_data/input.fasta.bz2 -o output_train_bz2 -m 2500 --ratio 0.05 -p 1')

assert os.path.exists('output_train_bz2/model.h5')

### Input .xz
os.system('SemiBin train --data test/train_data/data.csv --data-split test/train_data/data_split.csv -c test/train_data/cannot.txt --epoches 1 --batch-size 2048 --mode single -i test/train_data/input.fasta.xz -o output_train_xz -m 2500 --ratio 0.05 -p 1')

assert os.path.exists('output_train_xz/model.h5')

### train several samples
os.system('SemiBin train --data test/train_data/data.csv test/train_data/data.csv test/train_data/data.csv --data-split test/train_data/data_split.csv test/train_data/data_split.csv test/train_data/data_split.csv -c test/train_data/cannot.txt test/train_data/cannot.txt test/train_data/cannot.txt --epoches 1 --batch-size 2048 --mode several -i test/train_data/input.fasta.xz test/train_data/input.fasta.xz test/train_data/input.fasta.xz -o output_train_several_xz -m 2500 --ratio 0.05 -p 1')

assert os.path.exists('output_train_several_xz/model.h5')