from SemiBin.main import parse_args
from SemiBin.utils import validate_normalize_args
import logging

def test_parse_args():
    args = parse_args(
            ['single_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output']
            )
    validate_normalize_args(logging, args)
    assert args.training_type == 'semi'
    assert args.sequencing_type == 'short_read'

def test_parse_args_backcompat():
    for t in ['semi', 'self']:
        args = parse_args(
                ['single_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output',
                    '--training-type', t]
                )
        validate_normalize_args(logging, args)
        assert args.training_type == t
        assert args.sequencing_type == 'short_read'

        args2 = parse_args(
                ['single_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output',
                    f'--{t}-supervised']
                )
        validate_normalize_args(logging, args2)
        assert args2.training_type == t
        assert args2.sequencing_type == 'short_read'

    for mode_spec in ['--mode=several', '--train-from-many']:
        args = parse_args(
            ['train', mode_spec,
                 '--data', 'test/train_data/data.csv',
                 '--data-split', 'test/train_data/data_split.csv',
                 '-c', 'test/train_data/cannot.txt',
                 '-i', 'test/train_data/input.fasta',
                 '-o', 'test-outputs/output_train_fa',
                 '-p', '1'])
        validate_normalize_args(logging, args)
        assert args.mode == 'several'

    args = parse_args(
        ['train', '--mode=single',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-c', 'test/train_data/cannot.txt',
             '-i', 'test/train_data/input.fasta',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'])
    validate_normalize_args(logging, args)
    assert args.mode == 'single'

    args = parse_args(
        ['train',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-c', 'test/train_data/cannot.txt',
             '-i', 'test/train_data/input.fasta',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'])
    validate_normalize_args(logging, args)
    assert args.mode == 'single'
