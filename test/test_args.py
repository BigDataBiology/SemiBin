import pytest
from SemiBin.main import parse_args
from SemiBin.utils import validate_normalize_args
import logging

def test_parse_args():
    args = parse_args(
            ['single_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'],
            )
    validate_normalize_args(logging, args)
    assert args.training_type == 'self'
    assert args.sequencing_type == 'short_read'

def test_sequencing_type():
    args = parse_args(
            ['single_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '--sequencing-type', 'Long-Reads',
                '-o', 'output'],
            )
    validate_normalize_args(logging, args)
    assert args.training_type == 'self'
    assert args.sequencing_type == 'long_read'

def test_parse_args_backcompat():
    for t in ['semi', 'self']:
        args2 = parse_args(
                ['single_easy_bin',
                    '-i' ,'./test/single_sample_data/input.fasta',
                    '-b', './test/single_sample_data/input.sorted.bam',
                    '-o', 'output',
                    f'--{t}-supervised'],
                )
        validate_normalize_args(logging, args2)
        assert args2.training_type == t
        assert args2.sequencing_type == 'short_read'

    args = parse_args(
        ['train_semi', '--train-from-many',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-c', 'test/train_data/cannot.txt',
             '-i', 'test/train_data/input.fasta',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'],
        )
    validate_normalize_args(logging, args)
    assert args.mode == 'several'

    args = parse_args(
        ['train_semi',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-c', 'test/train_data/cannot.txt',
             '-i', 'test/train_data/input.fasta',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'])
    validate_normalize_args(logging, args)
    assert args.mode == 'single'

    args = parse_args(
        ['train_self',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-o', 'test-outputs/output_train_fa',
             '-p', '1'])
    validate_normalize_args(logging, args)
    assert args.mode == 'single'

    args = parse_args(
        ['train_self',
             '--data', 'test/train_data/data.csv',
             '--data-split', 'test/train_data/data_split.csv',
             '-o', 'test-outputs/output_train_fa',
             '--train-from-many',
             '-p', '1'])
    validate_normalize_args(logging, args)
    assert args.mode == 'several'


def test_single_easy_bin_abundance_args():
    args = parse_args(
            ['single_easy_bin',
                '-i', './test/coassembly_sample_data/input.fasta',
                '-a',
                './test/coassembly_sample_data/sample1.txt',
                './test/coassembly_sample_data/sample2.txt',
                './test/coassembly_sample_data/sample3.txt',
                './test/coassembly_sample_data/sample4.txt',
                './test/coassembly_sample_data/sample5.txt',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert args.training_type == 'self'
    assert args.sequencing_type == 'short_read'


def test_single_easy_bin_abundance_args_supervised():
    for t in ['semi', 'self']:
        args = parse_args(
                ['single_easy_bin',
                    f'--{t}-supervised',
                    '-i', './test/coassembly_sample_data/input.fasta',
                    '-a',
                    './test/coassembly_sample_data/sample1.txt',
                    './test/coassembly_sample_data/sample2.txt',
                    './test/coassembly_sample_data/sample3.txt',
                    './test/coassembly_sample_data/sample4.txt',
                    './test/coassembly_sample_data/sample5.txt',
                    '-o', 'output'])
        validate_normalize_args(logging, args)
        assert args.training_type == t


def test_quiet_before_or_after():
    args = parse_args([
            '--quiet',
            'single_easy_bin',
            '-i' ,'./test/single_sample_data/input.fasta',
            '-b', './test/single_sample_data/input.sorted.bam',
            '-o', 'output'])
    assert args.quiet
    args_after = parse_args([
            'single_easy_bin',
            '--quiet',
            '-i' ,'./test/single_sample_data/input.fasta',
            '-b', './test/single_sample_data/input.sorted.bam',
            '-o', 'output'])
    assert args_after.quiet
    assert args == args_after


def test_multi_easy_bin_args():
    args = parse_args(
            ['multi_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    assert args.cmd == 'multi_easy_bin'

def test_semibin2_args():
    args = parse_args(
            ['single_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert args.cmd == 'single_easy_bin'
    assert args.training_type == 'self'

    args = parse_args(
            ['single_easy_bin',
                '--self-supervised',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert args.cmd == 'single_easy_bin'
    assert args.training_type == 'self'

    args = parse_args(
            ['multi_easy_bin',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'],
            )
    validate_normalize_args(logging, args)
    assert args.cmd == 'multi_easy_bin'
    assert args.training_type == 'self'


    args = parse_args(
            ['multi_easy_bin',
                '--semi-supervised',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert args.cmd == 'multi_easy_bin'
    assert args.training_type == 'semi'

def test_write_prerecluster():
    args = parse_args(
            ['multi_easy_bin',
                '--semi-supervised',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert not bool(args.write_pre_reclustering_bins)

    args = parse_args(
            ['multi_easy_bin',
                '--semi-supervised',
                '--write-pre-reclustering-bins',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert args.write_pre_reclustering_bins

    args = parse_args(
            ['multi_easy_bin',
                '--semi-supervised',
                '--no-write-pre-reclustering-bins',
                '-i' ,'./test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'])
    validate_normalize_args(logging, args)
    assert not args.write_pre_reclustering_bins

    args = parse_args(
            ['multi_easy_bin',
                '--semi-supervised',
                '--write-pre-reclustering-bins',
                '--no-recluster',
                '-i', './test/single_sample_data/input.fasta',
                '-b', './test/single_sample_data/input.sorted.bam',
                '-o', 'output'],
            )
    with pytest.raises(SystemExit):
        validate_normalize_args(logging, args)


def test_bibtex():
    args = parse_args(['citation', '--bibtex'])
    assert args.cmd == 'citation'
    assert args.cite_format == 'bibtex'


def test_install_skills_args():
    args = parse_args(['install-skills'])
    assert args.cmd == 'install-skills'
    assert args.install_user is False
    assert args.skills_dir is None

    assert parse_args(['install-skills', '--user']).install_user is True
    assert parse_args(['install-skills', '--global']).install_user is True
    assert parse_args(['install-skills', '--skills-dir', '/tmp/x']).skills_dir == '/tmp/x'

    # the install_skills (underscore) alias normalizes to the canonical name
    assert parse_args(['install_skills']).cmd == 'install-skills'


def test_install_skills_installs(tmp_path):
    from SemiBin.main import install_skills
    dest = tmp_path / 'skills'
    install_skills(logging, to_user=False, skills_dir=str(dest))

    skill = dest / 'semibin' / 'SKILL.md'
    assert skill.exists()
    assert 'name: semibin' in skill.read_text()

    # Re-installing over an existing copy succeeds (over-writes in place)
    install_skills(logging, to_user=False, skills_dir=str(dest))
    assert skill.exists()
