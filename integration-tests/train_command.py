import os
import subprocess

### train from one sample
subprocess.check_call(
    ['SemiBin', 'train',
     '--data', 'test/train_data/data.csv',
     '--data-split', 'test/train_data/data_split.csv',
     '-c', 'test/train_data/cannot.txt',
     '--epoches', '1',
     '--batch-size', '2048',
     '-i', 'test/train_data/input.fasta',
     '-o', 'test-outputs/output_train_fa',
     '-m', '2500',
     '--ratio', '0.05',
     '-p', '1'])
assert os.path.exists('test-outputs/output_train_fa/model.h5')

subprocess.check_call(
    ['SemiBin', 'train_self',
     '--data', 'test/train_data/data.csv',
     '--data-split', 'test/train_data/data_split.csv',
     '--epoches', '1',
     '--batch-size', '2048',
     '-o', 'test-outputs/output_train_fa_self',
     '-p', '1'])
assert os.path.exists('test-outputs/output_train_fa_self/model.h5')

### train from several samples
subprocess.check_call(
        ['SemiBin', 'train',
         '--data',
             'test/train_data/data.csv',
             'test/train_data/data.csv',
             'test/train_data/data.csv',
         '--data-split',
             'test/train_data/data_split.csv',
             'test/train_data/data_split.csv',
             'test/train_data/data_split.csv',
         '-c',
             'test/train_data/cannot.txt',
             'test/train_data/cannot.txt',
             'test/train_data/cannot.txt',
         '-i',
             'test/train_data/input.fasta.gz',
             'test/train_data/input.fasta.xz',
             'test/train_data/input.fasta.bz2',
         '-o', 'test-outputs/output_train_several',
         '--epoches', '1',
         '--batch-size', '2048',
         '--train-from-many',
         '-m', '2500',
         '--ratio', '0.05',
         '-p', '1'])
assert os.path.exists('test-outputs/output_train_several/model.h5')

subprocess.check_call(
        ['SemiBin', 'train_self',
         '--data',
             'test/train_data/data.csv',
             'test/train_data/data.csv',
             'test/train_data/data.csv',
         '--data-split',
             'test/train_data/data_split.csv',
             'test/train_data/data_split.csv',
             'test/train_data/data_split.csv',
         '-o', 'test-outputs/output_train_several_self',
         '--epoches', '1',
         '--batch-size', '2048',
         '--train-from-many',
         '-p', '1'])
assert os.path.exists('test-outputs/output_train_several_self/model.h5')

