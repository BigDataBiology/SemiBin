from SemiBin.main import parse_args
from SemiBin.utils import validate_normalize_args
import logging

def test_parse_args():
    args = parse_args(
            ['single_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'],
            is_semibin2=False
            )
    validate_normalize_args(logging, args)
    assert args.training_type == 'semi'
    assert args.sequencing_type == 'short_read'

def test_sequencing_type():
    args = parse_args(
            ['single_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '--sequencing-type', 'Long-Reads',
                '-o', 'output'],
            is_semibin2=False
            )
    validate_normalize_args(logging, args)
    assert args.training_type == 'semi'
    assert args.sequencing_type == 'long_read'

def test_parse_args_backcompat():
    for t in ['semi', 'self']:
        args = parse_args(
                ['single_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output',
                    '--training-type', t],
                is_semibin2=False
                )
        validate_normalize_args(logging, args)
        assert args.training_type == t
        assert args.sequencing_type == 'short_read'

        args2 = parse_args(
                ['single_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output',
                    f'--{t}-supervised'],
                is_semibin2=False
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
                 '-p', '1'],
            is_semibin2=False,
            )
        validate_normalize_args(logging, args)
        assert args.mode == 'several'

    args = parse_args(
        ['train', '--mode=single',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-c', 'test/train_data/cannot.txt',
             '-i', 'test/train_data/input.fasta',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'],
        is_semibin2=False)
    validate_normalize_args(logging, args)
    assert args.mode == 'single'

    args = parse_args(
        ['train',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-c', 'test/train_data/cannot.txt',
             '-i', 'test/train_data/input.fasta',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'],
        is_semibin2=False)
    validate_normalize_args(logging, args)
    assert args.mode == 'single'

    args = parse_args(
        ['train_self',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'],
        is_semibin2=True)
    validate_normalize_args(logging, args)
    assert args.mode == 'single'

    args = parse_args(
        ['train_self',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-o', 'test-outputs/output_train_fa',
             '--train-from-many',
             '-p', '1'],
        is_semibin2=True)
    validate_normalize_args(logging, args)
    assert args.mode == 'several'


def test_quiet_before_or_after():
    args = parse_args([
            '--quiet',
            'single_easy_bin',
            '-i' ,'./test/single_sample_data/input.fasta',
            '-b', './test/single_sample_data/input.sorted.bam',
            '-o', 'output'],
            is_semibin2=False)
    assert args.quiet
    args_after = parse_args([
            'single_easy_bin',
            '--quiet',
            '-i' ,'./test/single_sample_data/input.fasta',
            '-b', './test/single_sample_data/input.sorted.bam',
            '-o', 'output'],
            is_semibin2=False)
    assert args_after.quiet
    assert args == args_after


def test_multi_easy_bin_args():
    args = parse_args(
            ['multi_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'],
            is_semibin2=False)
    assert args.cmd == 'multi_easy_bin'

def test_semibin2_args():
    for is_semibin2, def_mode in [(True, 'self'), (False, 'semi')]:
        args = parse_args(
                ['single_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert args.cmd == 'single_easy_bin'
        assert args.training_type == def_mode

        args = parse_args(
                ['single_easy_bin',
                    '--self-supervised',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert args.cmd == 'single_easy_bin'
        assert args.training_type == 'self'

        args = parse_args(
                ['multi_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert args.cmd == 'multi_easy_bin'
        assert args.training_type == def_mode


        args = parse_args(
                ['multi_easy_bin',
                    '--semi-supervised',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert args.cmd == 'multi_easy_bin'
        assert args.training_type == 'semi'

def test_write_prerecluster():
    for is_semibin2 in [True, False]:
        args = parse_args(
                ['multi_easy_bin',
                    '--semi-supervised',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert bool(args.write_pre_reclustering_bins) != is_semibin2

        args = parse_args(
                ['multi_easy_bin',
                    '--semi-supervised',
                    '--write-pre-reclustering-bins',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert args.write_pre_reclustering_bins

        args = parse_args(
                ['multi_easy_bin',
                    '--semi-supervised',
                    '--no-write-pre-reclustering-bins',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output'],
                is_semibin2=is_semibin2)
        validate_normalize_args(logging, args)
        assert not args.write_pre_reclustering_bins
