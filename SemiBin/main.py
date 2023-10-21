import argparse
import logging
import os
from multiprocessing.pool import Pool
import subprocess
from .atomicwrite import atomic_write
import shutil
import sys
from itertools import groupby
from .utils import validate_normalize_args, get_must_link_threshold, generate_cannot_link, \
    set_random_seed, process_fasta, split_data, get_model_path
from .generate_coverage import generate_cov, combine_cov
from .generate_kmer import generate_kmer_features_from_fasta
from .fasta import fasta_iter


def parse_args(args, is_semibin2):
    from .semibin_version import __version__
    # BooleanOptionalAction is available in Python 3.9; before that, we fall back on the default
    BooleanOptionalAction = getattr(argparse, 'BooleanOptionalAction', 'store_true')

    deprecated_if2_text = '[deprecated] ' if is_semibin2 else ''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description='Neural network-based binning of metagenomic contigs',
                                    epilog='For more information, see https://semibin.readthedocs.io/en/latest/subcommands/')

    parser.version = __version__

    parser.add_argument('-v',
                        '-V',
                        '--version',
                        action='version',
                        help='Print the version number')

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument('--verbose',
                   required=False,
                   help='Verbose output',
                   dest='verbose',
                   action='store_true', )

    verbosity.add_argument('--quiet', '-q',
                        required=False,
                        help='Quiet output',
                        dest='quiet',
                        action='store_true', )
    subparsers = parser.add_subparsers(title='SemiBin subcommands',
                                       dest='cmd',
                                       metavar='')

    single_easy_bin = subparsers.add_parser('single_easy_bin',
                                            help='Bin contigs (single or co-assembly) using one command.')

    multi_easy_bin = subparsers.add_parser('multi_easy_bin',
                                            help='Bin contigs (multi-sample mode) using one command.')


    generate_sequence_features_single = subparsers.add_parser('generate_sequence_features_single',
                                            help='Generate sequence features (kmer and abundance) as training data'
                                                  ' for (semi/self)-supervised deep learning model training (single or co-assembly mode).'
                                                  ' This will produce the data.csv and data_split.csv files.'
                                                  )


    generate_sequence_features_multi = subparsers.add_parser('generate_sequence_features_multi', aliases=['generate_sequence_features_multi'],
                                            help='Generate sequence features (kmer and abundance) as training data'
                                                  ' for (semi/self)-supervised deep learning model training (multi-sample mode).'
                                                  ' This will produce the data.csv and data_split.csv files.'
                                                  )


    check_install = subparsers.add_parser('check_install', help = 'Check whether required dependencies are present.')

    check_install.add_argument('--allow-missing-mmseqs2',
            required=False,
            help='Do not fail is MMSeqs2 is not found. MMSeqs2 is required for semi-supervised learning, but not self-supervised learning.',
            dest='allow_missing_mmseqs2',
            action='store_true', )

    concatenate_fasta = subparsers.add_parser('concatenate_fasta', help = 'concatenate fasta files for multi-sample binning')

    concatenate_fasta.add_argument('-i', '--input-fasta',
                   required=True,
                   nargs='*',
                   help='Path to the input fasta file.',
                   dest='contig_fasta',
                   default=None, )

    concatenate_fasta.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None
                        )

    concatenate_fasta.add_argument('-m', '--min-len',
                        required=False,
                        type=int,
                        help='Discard sequences below this length (default:0)',
                        default=0,
                        dest='min_length')

    generate_sequence_features_single.add_argument('--kmer',
                                                   required=False,
                                                   help='Just output data.csv with k-mer features.',
                                                   dest='kmer',
                                                   action='store_true',)


    training_self = subparsers.add_parser('train_self',
                                          help = 'Train the model with self-supervised learning')

    binning = subparsers.add_parser('bin', aliases=['bin_short'],
                                    help='Group the contigs into bins.')

    binning_long = subparsers.add_parser('bin_long',
                                    help='Group the contigs from long reads into bins.')


    # Add the deprecated arguments last so they show up at the bottom of the help text
    train_semi = subparsers.add_parser(('train_semi' if is_semibin2 else 'train'),
                                    help=deprecated_if2_text + 'Train the model.')

    generate_cannot_links = subparsers.add_parser('generate_cannot_links', aliases=['predict_taxonomy',],
                                             help=deprecated_if2_text + 'Run the contig annotation using mmseqs '
                                                  'with GTDB reference genome and generate '
                                                  'cannot-link file used in the semi-supervised deep learning model training. '
                                                  'This will download the GTDB database if not downloaded before.')
    download_GTDB = subparsers.add_parser('download_GTDB',
                help=deprecated_if2_text + 'Download GTDB reference genomes.')

    download_GTDB.add_argument('-f', '--force',
                            required=False,
                            help='Redownload GTDB even if files are found',
                            dest='force',
                            action='store_true',
                            default=None)

    training_mandatory = train_semi.add_argument_group('Mandatory arguments')
    training_mandatory.add_argument('-i', '--input-fasta',
                   required=True,
                   nargs='*',
                   help='Path to the input fasta file.',
                   dest='contig_fasta',
                   default=None, )
    training_mandatory.add_argument('-c', '--cannot-link',
                         required=True,
                          nargs='*',
                         help='Path to the input cannot link file. '
                         'The file format: `contig_1,contig_2` '
                         '(one row for each cannot link constraint).',
                         dest='cannot_link',
                         default=None)

    training_self_mandatory = training_self.add_argument_group('Mandatory arguments')
    for p in [training_mandatory, training_self_mandatory]:
        p.add_argument('--data',
                             required=True,
                             nargs='*',
                             help='Path to the input data.csv file.',
                             dest='data',
                             default=None,
                             )
        p.add_argument('--data-split',
                             required=True,
                             nargs='*',
                             help='Path to the input data_split.csv file.',
                             dest='data_split',
                             default=None,
                             )
        p.add_argument('-o', '--output',
                              required=True,
                              help='Output directory (will be created if non-existent)',
                              dest='output',
                              default=None,
                              )

    for p in [train_semi, training_self]:
        p.add_argument('--batch-size',
                              required=False,
                              type=int,
                              help='Batch size used in the training process (Default: 2048).',
                              dest='batchsize',
                              default=2048, )

        if not is_semibin2:
            p.add_argument('--mode',
                              required=False,
                              type=str,
                              help='[DEPRECATED: use --train-from-many]. [single/several] Train models from one (single) or more samples (several). '
                                   'In `several` mode, you must provide data, data_split, cannot, and fasta files for corresponding samples in the same order. '
                                   'Note: You can only use `several` mode when performing single-sample binning. Training from several samples with multi-sample binning is not supported.',
                              dest='mode',
                              default='single')

        p.add_argument('--train-from-many',
                           required=False,
                           help='Train the model with several samples.\n'
                                   'You must provide data, data_split, cannot, and fasta files for corresponding samples in the same order. '
                                   'Note: You can only use `--train-from-many` mode when performing single-sample binning. Training from many samples with multi-sample binning is not supported.',
                           dest='train_from_many',
                           action=BooleanOptionalAction)



    train_semi.add_argument('--epochs', '--epoches', # epoches is kept for backwards compatibilty
                       required=False,
                       type=int,
                       help='Number of epochs used in the training process (Default: 20).',
                       dest='epoches',
                       default=20)


    training_self.add_argument('--epochs', '--epoches', # epoches is kept for backwards compatibilty
                       required=False,
                       type=int,
                       help='Number of epochs used in the training process (Default: 15).',
                       dest='epoches',
                       default=15)



    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links,
                generate_sequence_features_single, generate_sequence_features_multi,
                binning, binning_long]:
        m = p.add_argument_group('Mandatory arguments')

        m.add_argument('-i', '--input-fasta',
                                required=True,
                                help='Path to the input fasta file.',
                                dest='contig_fasta',
                                default=None,)
        m.add_argument('-o', '--output',
                            required=True,
                            help='Output directory (will be created if non-existent)',
                            dest='output',
                            default=None,
                            )
        if p is binning or p is binning_long:
            m.add_argument('--data',
                             required=True,
                             help='Path to the input data.csv file.',
                             dest='data',
                             default=None,)
        if p in [multi_easy_bin, generate_sequence_features_multi]:
            m.add_argument('-b', '--input-bam',
                              required=True,
                              nargs='*',
                              help='Path to the input BAM(.bam)/CRAM(.cram) file. '
                                   'If using multiple samples, you can input multiple files.',
                              dest='bams',
                              default=None,
                              )
        p.add_argument('--write-pre-reclustering-bins',
                required=False,
                help='Write pre-reclustering bins to disk.',
                dest='write_pre_reclustering_bins',
                action=BooleanOptionalAction)
        if not hasattr(argparse, 'BooleanOptionalAction'):
            p.add_argument('--no-write-pre-reclustering-bins',
                    required=False,
                    help='Do not write pre-reclustering bins to disk.',
                    dest='no_write_pre_reclustering_bins',
                    action='store_true')

        p.add_argument('--tag-output',
                required=False,
                type=str,
                dest='output_tag',
                default=('SemiBin' if is_semibin2 else None),
                help='Tag to add to output file names')

    for p in [single_easy_bin,
                multi_easy_bin,
                generate_cannot_links,
                generate_sequence_features_single,
                generate_sequence_features_multi,
                binning, binning_long,
                concatenate_fasta]:
        p.add_argument('--compression',
                required=False,
                type=str,
                help=('Compression type for the output files (accepted values: ' +
                    ('none [default]/gz/xz/bz2).' if not is_semibin2 else
                    ' gz [default]/xz/bz2/none).')),
                dest='output_compression',
                default=('gz' if is_semibin2 else 'none'))

    for p in [binning, binning_long]:
        p.add_argument('--model',
                             required=False,
                             type=str,
                             dest='model_path',
                             default=None,
                             help='Path to the trained deep learning model.')


    for p in [single_easy_bin, multi_easy_bin, train_semi, binning, binning_long]:
        p.add_argument('--orf-finder',
                       required=False,
                       type=str,
                       help='ORF finder used to estimate the number of bins (prodigal/fraggenescan)',
                       dest='orf_finder',
                       default=('fast-naive' if is_semibin2 else 'prodigal'))
        p.add_argument('--prodigal-output-faa',
                       required=False,
                       type=str,
                       help='Bypasses ORF calling and uses the provided .faa file instead (must be in same format as prodigal output).',
                       dest='prodigal_output_faa')

    for p in [single_easy_bin, binning, binning_long]:
        p.add_argument('--depth-metabat2',
                       required=False,
                       type=str,
                       help='depth file generated by metabat2 (only used with single-sample binning)',
                       dest='depth_metabat2',
                       default=None,
                       )

    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links, train_semi, binning, generate_sequence_features_single, generate_sequence_features_multi, training_self, binning_long]:
        p.add_argument('--tmpdir',
                       required=False,
                       type=str,
                       help='option to set temporary directory',
                       dest='tmp_output',
                       default=None,
                       )


    for p in [train_semi, generate_cannot_links, binning, single_easy_bin, multi_easy_bin, generate_sequence_features_single, generate_sequence_features_multi, training_self, binning_long]:
        p.add_argument('-p', '--processes', '-t', '--threads',
                   required=False,
                   type=int,
                   help='Number of CPUs used (pass the value 0 to use all CPUs, default: 0)',
                   dest='num_process',
                   default=0,
                   metavar=''
                   )

    for p in [single_easy_bin, binning, binning_long]:
        p.add_argument('--environment', '--habitat', '--biome',
                       required=False,
                       help='Environment for the built-in model (available choices: human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global).',
                       dest='environment',
                       default=None,
                       )

    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links, generate_sequence_features_single, generate_sequence_features_multi, binning, train_semi, binning_long]:
        p.add_argument('-m', '--min-len',
                       required=False,
                       type=int,
                       help='Minimal length for contigs in binning. '
                            'If you use SemiBin with multi steps and you use this parameter, please use this parameter consistently with all subcommands. '
                            '(Default: SemiBin chooses 1000bp or 2500bp according the ratio of the number of base pairs of contigs between 1000-2500bp).',
                       dest='min_len',
                       default=None,
                       )
        p.add_argument('--ratio',
                       required=False,
                       type=float,
                       help='If the ratio of the number of base pairs of contigs between 1000-2500 bp smaller than this value, the minimal length will be set as 1000bp, otherwise 2500bp. '
                       'Note that setting `--min-length/-m` overrules this parameter. '
                       'If you use SemiBin with multi steps and you use this parameter, please use this parameter consistently with all subcommands. '
                       '(Default: 0.05)',
                       dest='ratio',
                       default=0.05)

    for p in [single_easy_bin, generate_sequence_features_single]:
        p.add_argument('-b', '--input-bam',
                              required=False,
                              nargs='*',
                              help='Path to the input BAM(.bam)/CRAM(.cram) file. '
                                   'If using multiple samples, you can input multiple files.'
                                    'If just need k-mer features, bam file is not needed.',
                              dest='bams',
                              default=None,
                              )

    for p in [single_easy_bin, multi_easy_bin,
                    generate_sequence_features_single, generate_sequence_features_multi,
                    generate_cannot_links, check_install, concatenate_fasta,
                    train_semi, binning, training_self, binning_long]:
        verbosity = p.add_mutually_exclusive_group()
        # Using verbose1/quiet1 is a hack for the fact that it is hard to make
        # argparse accept options both in the global scope and in the
        # subcommand scope
        verbosity.add_argument('--verbose',
                        required=False,
                        help='Verbose output',
                        dest='verbose1',
                        action='store_true', )
        verbosity.add_argument('--quiet', '-q',
                        required=False,
                        help='Quiet output',
                        dest='quiet1',
                        action='store_true', )
    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links, download_GTDB]:
        p.add_argument('-r', '--reference-db-data-dir', '--reference-db',
                            required=False,
                            help='GTDB reference storage path. (Default: $HOME/.cache/SemiBin/mmseqs2-GTDB).'
                            'If not set --reference-db and SemiBin cannot find GTDB in $HOME/.cache/SemiBin/mmseqs2-GTDB, SemiBin will download GTDB (Note that >100GB of disk space are required).',
                            dest='GTDB_reference',
                            metavar='',
                            default=None)

    for p in [single_easy_bin, generate_cannot_links, multi_easy_bin]:
        p.add_argument('--cannot-name',
                            required=False,
                            help='Name for the cannot-link file (default: cannot).',
                            dest='cannot_name',
                            default='cannot',
                            metavar=''
                            )
        p.add_argument('--taxonomy-annotation-table',
                            required=False,
                            nargs='*',
                            help='Pre-computed mmseqs2 format taxonomy TSV file to bypass mmseqs2 GTDB annotation [advanced]. '
                                 'When running with multi-sample binning, please make sure that the order of the taxonomy TSV file and the contig file (used for the combined fasta) is same.',
                            dest='taxonomy_results_fname',
                            metavar='TAXONOMY_TSV')

    binning_long.add_argument('--minfasta-kbs',
                            required=False,
                            type=int,
                            help='minimum bin size in Kbps (Default: 200).',
                            dest='minfasta_kb',
                            default=200,
                            metavar='')

    for p in [single_easy_bin, multi_easy_bin]:
        p.add_argument('--epochs', '--epoches', # epoches is kept for backwards compatibilty
                       required=False,
                       type=int,
                       help='Number of epochs used in the training process (Default: 15).',
                       dest='epoches',
                       default=15)

        p.add_argument('--batch-size',
                       required=False,
                       type=int,
                       help='Batch size used in the training process (Default: 2048).',
                       dest='batchsize',
                       default=2048)


    for p in [single_easy_bin, multi_easy_bin, binning]:
        p.add_argument('--minfasta-kbs',
                            required=False,
                            type=int,
                            help='minimum bin size in Kbps (Default: 200).',
                            dest='minfasta_kb',
                            default=200,
                            metavar='')

        g = p.add_argument_group('Binning options (advanced use)')
        g.add_argument('--max-edges',
                          required=False,
                          type=int,
                          help='The maximum number of edges that can be connected to one contig (Default: 200).',
                          dest='max_edges',
                          default=200)

        g.add_argument('--max-node',
                          required=False,
                          type=float,
                          dest='max_node',
                          default=1,
                          help='Fraction of contigs that considered to be binned (should be between 0 and 1; default: 1).')

        g.add_argument('--no-recluster',
                           required=False,
                           help='Do not recluster bins.',
                           dest='no_recluster',
                           action='store_true', )

        if not is_semibin2:
            g.add_argument('--recluster',
                           required=False,
                           help='[Deprecated] Does nothing (current default is to perform clustering)',
                           dest='recluster',
                           action='store_true', )


    for p in [multi_easy_bin, generate_sequence_features_multi, concatenate_fasta]:
        p.add_argument('-s', '--separator',
                           required=False,
                           type=str,
                           help='Used when multiple samples binning to separate sample name and contig name (Default is :).',
                           dest='separator',
                           default=':',
                           metavar='')

    for p in [train_semi, binning, single_easy_bin, multi_easy_bin, training_self, binning_long]:
        p.add_argument('--random-seed',
                       required=False,
                       type=int,
                       help='Random seed. Set it to a fixed value to reproduce results across runs. The default is that the seed is set by the system and .',
                       dest='random_seed',
                       default=None,
                       )

    for p in [generate_cannot_links, generate_sequence_features_single, generate_sequence_features_multi, single_easy_bin, multi_easy_bin]:
        p.add_argument('--ml-threshold',
                       required=False,
                       type=int,
                       help='Length threshold for generating must-link constraints. (By default, the threshold is calculated from the contig, and the default minimum value is 4,000 bp)',
                       dest='ml_threshold',
                       default=None)

    for p in [single_easy_bin, multi_easy_bin, train_semi, binning, training_self, binning_long]:
        p.add_argument('--engine',
                       required=False,
                       type=str,
                       help='device used to train the model (auto/gpu/cpu, auto means if SemiBin detects the gpu, SemiBin will use GPU)',
                       dest='engine',
                       default='auto')

    for p in [single_easy_bin, multi_easy_bin]:
        p.add_argument('--semi-supervised',
                           required=False,
                           help='Train the model with semi-supervised learning.',
                           dest='semi_supervised',
                           action='store_true', )

        p.add_argument('--self-supervised',
                           required=False,
                           help='Train the model with self-supervised learning.',
                           dest='self_supervised',
                           action='store_true', )

        if not is_semibin2:
            p.add_argument('--training-type',
                       required=False,
                       type=str,
                       help='Training algorithm used to train the model (semi [default]/self)\n'
                            'DEPRECATED: use --self-supervised/--semi-supervised',
                       dest='training_type')

        p.add_argument('--sequencing-type',
               required=False,
               type=str,
               help='sequencing type in [short_read/long_read], Default: short_read.',
               dest='sequencing_type',
               default='short_read',)


    if not args:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args(args)
    args.is_semibin2 = is_semibin2
    if hasattr(args, 'no_recluster'):
        args.recluster = not args.no_recluster

    if hasattr(args, 'write_pre_reclustering_bins'):
        # backwards compat for Python 3.8 and earlier
        if not hasattr(argparse, 'BooleanOptionalAction'):
            if args.no_write_pre_reclustering_bins:
                args.write_pre_reclustering_bins = False
            elif not args.no_write_pre_reclustering_bins and not is_semibin2:
                args.write_pre_reclustering_bins = True
        if args.write_pre_reclustering_bins is None:
            args.write_pre_reclustering_bins = not is_semibin2

    # Keep the verbose1/quiet1 hack contained in this function
    for hacked in ['verbose', 'quiet']:
        if hasattr(args, f'{hacked}1'):
            if getattr(args, f'{hacked}1'):
                setattr(args, hacked, getattr(args, f'{hacked}1'))
            delattr(args, f'{hacked}1')

    if args.cmd == 'bin_short':
        args.cmd = 'bin'
    return args



def check_install(verbose, orf_finder=None, allow_missing_mmseqs2=False):
    '''Executes check_install subcommand, which checks for dependencies

    Parameters
    ----------
    verbose : boolean. If true, then prints out all the paths
    orf_finder : str, optional
    allow_missing_mmseqs2 : boolean, optional
        If true, then checks for mmseqs2
    '''
    from shutil import which
    dependencies = ['bedtools', 'hmmsearch', 'mmseqs', 'FragGeneScan', 'prodigal']
    has_fgs = False
    missing_deps = False
    if verbose:
        print("Looking for dependencies...")
    for dep in dependencies:
        p = which(dep)
        if not p:
            if dep == 'mmseqs':
                if not allow_missing_mmseqs2:
                    sys.stderr.write(
                        f"Error: {dep} does not seem to be installed! This is necessary for semi-supervised learning\n")
                    missing_deps = True
                elif verbose:
                    print(f'\t{dep} not found. Semi-supervised training will not be possible')
            elif dep == 'prodigal':
                if orf_finder == 'fast-naive':
                    pass
                elif not has_fgs:
                    if orf_finder != 'fast-naive':
                        sys.stderr.write(
                                'Error: neither prodigal nor FragGeneScan appear to be available!\n'
                                'You can use --orf-finder=fast-naive to use the builtin simple ORF finder')
                        missing_deps = True
                else:
                    if verbose or orf_finder == 'prodigal':
                        sys.stderr.write(
                            'Warning: prodigal does not appear to be available (although FragGeneScan is). You must use the `--orf-finder=fast-naive` or `--orf-finder=fraggenescan` options.\n')
                    missing_deps = True
            elif dep == 'FragGeneScan':
                pass
            else:
                sys.stderr.write(
                    f"Error: {dep} does not seem to be installed!\n")
                missing_deps = True
        else:
            if dep == 'FragGeneScan':
                has_fgs = True
            if verbose:
                print(f'\t{dep:16}: {p}')
    if missing_deps:
        print('Missing dependencies')
        sys.exit(1)
    else:
        if verbose:
            print('Installation OK')


def predict_taxonomy(logger, contig_fasta, cannot_name,
                     taxonomy_results_fname, GTDB_reference, num_process,
                     binned_length, must_link_threshold, output):
    """
    Predict taxonomy using mmseqs and generate cannot-link file

    contig_fasta: contig file used
    GTDB_reference: GTDB path
    cannot_name: name for cannot-link constraint
    binned_short: threshold for contigs used in binning
    must_link_threshold: threshold of contigs for must-link constraints
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tdir:
        filtered_fasta = os.path.join(tdir, 'SemiBinFiltered.fa')
        namelist = []
        num_must_link = 0
        with open(filtered_fasta, 'wt') as out:
            for h,seq in fasta_iter(contig_fasta):
                if len(seq) >= binned_length:
                    out.write(f'>{h}\n{seq}\n')
                    namelist.append(h)
                    if len(seq) >= must_link_threshold:
                        num_must_link += 1
        if taxonomy_results_fname is None:
            from .gtdb import find_or_download_gtdb
            GTDB_reference = find_or_download_gtdb(logger, GTDB_reference, force=False)
            try:
                subprocess.check_call(
                    ['mmseqs',
                     'createdb',
                     filtered_fasta,
                     os.path.join(tdir, 'contig_DB')],
                    stdout=None,
                )
            except:
                sys.stderr.write(
                    f"Error: Running mmseqs createdb fail\n")
                sys.exit(1)
            if os.path.exists(os.path.join(output, 'mmseqs_contig_annotation')):
                shutil.rmtree(os.path.join(output, 'mmseqs_contig_annotation'))
            os.makedirs(os.path.join(output, 'mmseqs_contig_annotation'))
            try:
                subprocess.run(
                    ['mmseqs',
                     'taxonomy',
                     os.path.join(tdir, 'contig_DB'),
                     GTDB_reference,
                     os.path.join(output, 'mmseqs_contig_annotation/mmseqs_contig_annotation'),
                     tdir,
                     '--tax-lineage', '1',
                     '--threads', str(num_process),
                     ],
                    check=True,
                    stdout=None,
                )
            except:
                sys.stderr.write(
                    f"Error: Running mmseqs taxonomy fail\n")
                sys.exit(1)
            taxonomy_results_fname = os.path.join(output,
                                        'mmseqs_contig_annotation',
                                        'taxonomyResult.tsv')
            try:
                subprocess.check_call(
                    ['mmseqs',
                     'createtsv',
                     os.path.join(tdir, 'contig_DB'),
                     os.path.join(output, 'mmseqs_contig_annotation/mmseqs_contig_annotation'),
                     taxonomy_results_fname,
                     ],
                    stdout=None,
                )
            except:
                sys.stderr.write(
                    f"Error: Running mmseqs createtsv fail\n")
                sys.exit(1)

    os.makedirs(os.path.join(output, 'cannot'), exist_ok=True)
    generate_cannot_link(taxonomy_results_fname,
        namelist, num_must_link, os.path.join(output, 'cannot'), cannot_name)

def generate_sequence_features_single(logger, contig_fasta,
                         bams, binned_length,
                         must_link_threshold, num_process, output, only_kmer=False):
    """
    Generate data.csv and data_split.csv for training and clustering of single-sample and co-assembly binning mode.
    data.csv has the features(kmer and abundance) for original contigs.
    data_split.csv has the features(kmer and abundace) for contigs that are breaked up as must-link pair.

    """
    import pandas as pd

    if bams is None and not only_kmer:
        sys.stderr.write(
            f"Error: You need to specify input BAM files if you want to calculate coverage features.\n")
        sys.exit(1)

    if bams is not None and only_kmer:
        logger.info('We will only calculate k-mer features.')

    if not only_kmer:
        n_sample = len(bams)
        is_combined = n_sample >= 5
        bam_list = bams

        logger.info('Calculating coverage for every sample.')

        with Pool(num_process if num_process != 0 else None) as pool:
            results = [
                pool.apply_async(
                    generate_cov,
                    args=(
                        bam_file,
                        bam_index,
                        output,
                        must_link_threshold,
                        is_combined,
                        binned_length,
                        logger,
                        None
                    ))
                for bam_index, bam_file in enumerate(bams)]
            for r in results:
                s = r.get()
                logger.info(f'Processed: {s}')

        for bam_index, bam_file in enumerate(bams):
            if not os.path.exists(os.path.join(output, '{}_data_cov.csv'.format(
                    os.path.split(bam_file)[-1] + '_{}'.format(bam_index)))):
                sys.stderr.write(
                    f"Error: Generating coverage file fail\n")
                sys.exit(1)
            if is_combined:
                if not os.path.exists(os.path.join(output, '{}_data_split_cov.csv'.format(
                        os.path.split(bam_file)[-1] + '_{}'.format(bam_index)))):
                    sys.stderr.write(
                        f"Error: Generating coverage file fail\n")
                    sys.exit(1)

        logger.debug('Start generating kmer features from fasta file.')
        kmer_whole = generate_kmer_features_from_fasta(
            contig_fasta, binned_length, 4)
        kmer_split = generate_kmer_features_from_fasta(
            contig_fasta, 1000, 4, split=True, split_threshold=must_link_threshold)

        data = kmer_whole
        data_split = kmer_split
        data.index = data.index.astype(str)

        if is_combined:
            data_cov, data_split_cov = combine_cov(output, bam_list, is_combined)
            data_split = pd.merge(data_split, data_split_cov, how='inner', on=None,
                                      left_index=True, right_index=True, sort=False, copy=True)
        else:
            data_cov = combine_cov(output, bam_list, is_combined)
            data_split = kmer_split

        data = pd.merge(data, data_cov, how='inner', on=None,
                                      left_index=True, right_index=True, sort=False, copy=True)


        with atomic_write(os.path.join(output, 'data.csv'), overwrite=True) as ofile:
            data.to_csv(ofile)

        with atomic_write(os.path.join(output, 'data_split.csv'), overwrite=True) as ofile:
            data_split.to_csv(ofile)

    else:
        logger.info('Only generating kmer features from fasta file.')
        kmer_whole = generate_kmer_features_from_fasta(
            contig_fasta, binned_length, 4)
        with atomic_write(os.path.join(output, 'data.csv'), overwrite=True) as ofile:
            kmer_whole.to_csv(ofile)


def generate_sequence_features_multi(logger, contig_fasta,
                        bams, num_process,
                        separator, ratio, min_length, ml_threshold, output):
    """
    Generate data.csv and data_split.csv for every sample of multi-sample binning mode.
    data.csv has the features(kmer and abundance) for original contigs.
    data_split.csv has the features(kmer and abundace) for contigs that are breaked up as must-link pair.
    """
    import pandas as pd
    n_sample = len(bams)
    is_combined = n_sample >= 5
    bam_list = bams

    # Gererate contig file for every sample
    sample_list = []
    contig_length_list = []

    os.makedirs(os.path.join(output, 'samples'), exist_ok=True)

    def fasta_sample_iter(fn):
        for h,seq in fasta_iter(fn):
            if separator not in h:
                raise ValueError(
                    f"Expected contigs to contain separator character ({separator}), found {h}")
            sample_name, contig_name = h.split(separator, 1)
            yield sample_name, contig_name, seq

    for sample_name, contigs in groupby(fasta_sample_iter(contig_fasta), lambda sn_cn_seq : sn_cn_seq[0]):
        with open(os.path.join( output, 'samples', '{}.fa'.format(sample_name)), 'wt') as out:
            for _, contig_name, seq in contigs:
                out.write(f'>{contig_name}\n{seq}\n')
                contig_length_list.append(len(seq))
        sample_list.append(sample_name)

    must_link_threshold = get_must_link_threshold(contig_length_list) if ml_threshold is None else ml_threshold

    logger.info('Calculating coverage for every sample.')

    binning_threshold = {}
    for sample in sample_list:
        if min_length is not None:
            binning_threshold[sample] = min_length
        else:
            binned_short ,_ ,_ = process_fasta(os.path.join(output, 'samples/{}.fa'.format(sample)), ratio)
            binning_threshold[sample] = 1000 if binned_short else 2500

    with Pool(num_process if num_process != 0 else None) as pool:
        results = [
            pool.apply_async(
                        generate_cov,
                         args=(
                             bam_file,
                             bam_index,
                             os.path.join(output, 'samples'),
                             must_link_threshold,
                             is_combined,
                             binning_threshold,
                             logger,
                             separator,
                         ))
            for bam_index, bam_file in enumerate(bams)]
        for r in results:
            s = r.get()
            logger.info(f'Processed: {s}')

    for bam_index, bam_file in enumerate(bam_list):
        if not os.path.exists(os.path.join(os.path.join(output, 'samples'), '{}_data_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index)))):
            sys.stderr.write(
                f"Error: Generating coverage file fail\n")
            sys.exit(1)
        if is_combined:
            if not os.path.exists(os.path.join(os.path.join(output, 'samples'), '{}_data_split_cov.csv'.format(
                    os.path.split(bam_file)[-1] + '_{}'.format(bam_index)))):
                sys.stderr.write(
                    f"Error: Generating coverage file fail\n")
                sys.exit(1)

    # Generate cov features for every sample
    if is_combined:
        data_cov, data_split_cov = combine_cov(os.path.join(output, 'samples'), bam_list, is_combined)
        data_split_cov = data_split_cov.reset_index()
        columns_list = list(data_split_cov.columns)
        columns_list[0] = 'contig_name'
        data_split_cov.columns = columns_list
    else:
        data_cov = combine_cov(os.path.join(output, 'samples'), bam_list, is_combined)

    data_cov = data_cov.reset_index()
    columns_list = list(data_cov.columns)
    columns_list[0] = 'contig_name'
    data_cov.columns = columns_list

    for sample in sample_list:
        output_path = os.path.join(output, 'samples', sample)
        os.makedirs(output_path, exist_ok=True)

        part_data = split_data(data_cov, sample, separator, is_combined)
        part_data.to_csv(os.path.join(output_path, 'data_cov.csv'))

        if is_combined:
            part_data = split_data(data_split_cov, sample, separator, is_combined)
            part_data.to_csv(os.path.join(
                output_path, 'data_split_cov.csv'))

        sample_contig_fasta = os.path.join(
            output, 'samples/{}.fa'.format(sample))
        binned_length =  binning_threshold[sample]
        kmer_whole = generate_kmer_features_from_fasta(
            sample_contig_fasta, binned_length, 4)
        kmer_split = generate_kmer_features_from_fasta(
            sample_contig_fasta, 1000, 4, split=True, split_threshold=must_link_threshold)


        sample_cov = pd.read_csv(os.path.join(output_path, 'data_cov.csv'),index_col=0)
        kmer_whole.index = kmer_whole.index.astype(str)
        sample_cov.index = sample_cov.index.astype(str)
        data = pd.merge(kmer_whole, sample_cov, how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)
        if is_combined:
            sample_cov_split = pd.read_csv(os.path.join(
                output_path, 'data_split_cov.csv'), index_col=0)
            data_split = pd.merge(kmer_split, sample_cov_split, how='inner', on=None,
                                  left_index=True, right_index=True, sort=False, copy=True)
        else:
            data_split = kmer_split

        with atomic_write(os.path.join(output_path, 'data.csv'), overwrite=True) as ofile:
            data.to_csv(ofile)

        with atomic_write(os.path.join(output_path, 'data_split.csv'), overwrite=True) as ofile:
            data_split.to_csv(ofile)

    return sample_list


def training(logger, contig_fasta, num_process,
             data, data_split, cannot_link, batchsize,
             epoches,  output, device, ratio, min_length, mode, *,
             orf_finder=None, prodigal_output_faa=None, training_mode='semi'):
    """
    Training the model

    model: [single/several]
    """
    from .semi_supervised_model import train
    from .self_supervised_model import train_self
    import pandas as pd
    binned_lengths = []

    if mode == 'single':
        logger.info('Start training from one sample.')
        data_ = pd.read_csv(data[0], index_col=0)
        col_name = data_.columns.tolist()[-1].split('_')[-1]
        is_combined = False if col_name == 'var' else True

    else:
        logger.info('Start training from multiple samples.')
        is_combined = False

    if training_mode == 'semi':
        contig_fasta_unzip = []
        for fasta_index,temp_fasta in enumerate(contig_fasta):
            if temp_fasta.endswith('.gz') or temp_fasta.endswith('.bz2') or temp_fasta.endswith('.xz'):
                temp_fasta_unzip = os.path.join(output, 'unzip_contig_{}.fa'.format(fasta_index))
                with open(temp_fasta_unzip, 'wt') as out:
                    for h,seq in fasta_iter(temp_fasta):
                            out.write(f'>{h}\n{seq}\n')
                contig_fasta_unzip.append(temp_fasta_unzip)
            else:
                contig_fasta_unzip.append(temp_fasta)

        if mode == 'single':
            if min_length is None:
                binned_short, _, _ = process_fasta(contig_fasta[0], ratio)
                binned_lengths.append(1000) if binned_short else binned_lengths.append(2500)
            else:
                binned_lengths.append(min_length)
        else:
            for contig_index in contig_fasta:
                if min_length is None:
                    binned_short, _, _ = process_fasta(contig_index, ratio)
                    binned_lengths.append(1000) if binned_short else binned_lengths.append(2500)
                else:
                    binned_lengths.append(min_length)


        model = train(
            output,
            contig_fasta_unzip,
            binned_lengths,
            logger,
            data,
            data_split,
            cannot_link,
            is_combined,
            batchsize,
            epoches,
            device,
            num_process,
            mode = mode,
            prodigal_output_faa=prodigal_output_faa,
            orf_finder=orf_finder)
    else:
        model = train_self(output,
                           logger,
                           data,
                           data_split,
                           is_combined,
                           batchsize,
                           epoches,
                           device,
                           num_process,
                           mode)


def binning_preprocess(data, depth_metabat2, model_path, environment, device):
    import pandas as pd
    import torch
    data = pd.read_csv(data, index_col=0)
    data.index = data.index.astype(str)

    if depth_metabat2 is None:
        col_name = data.columns.tolist()[-1].split('_')[-1]
        is_combined = col_name != 'var'
        n_sample = (data.shape[1] - 136) // 2 if not is_combined else (data.shape[1] - 136)
    else:
        is_combined = False
        n_sample = 1
        depth_metabat2 = pd.read_csv(depth_metabat2, sep='\t')
        depth_metabat2.set_index(['contigName'], inplace=True)
        depth_metabat2.drop(['contigLen', 'totalAvgDepth'], axis=1, inplace=True)
        depth_metabat2.index.name = None
        data = pd.merge(data, depth_metabat2, how='inner', on=None,
                 left_index=True, right_index=True, sort=False, copy=True)
        if data.shape[1] != 138:
            sys.stderr.write(
                f"Error: Depth file from Metabat2 can only be used in single-sample binning!\n")
            sys.exit(1)

    model_path = model_path if environment is None else get_model_path(environment)
    if environment is not None:
        if data.shape[1] != 138:
            sys.stderr.write(f"Error: provided pretrained model only used in single-sample binning!\n")
            sys.exit(1)

    if device == torch.device('cpu'):
        model = torch.load(model_path, map_location=torch.device('cpu'))
    else:
        model = torch.load(model_path).to(device)

    return is_combined, n_sample, data, model

def binning_long(logger, data, minfasta, binned_length, contig_dict,
        model_path, output, device, environment, *, args):
    from .long_read_cluster import cluster_long_read
    logger.info('Start binning.')
    is_combined, n_sample, data, model = binning_preprocess(data, getattr(args, 'depth_metabat2', None), model_path, environment, device)
    cluster_long_read(model,
                      data,
                      device,
                      is_combined,
                      logger,
                      n_sample,
                      output,
                      contig_dict,
                      binned_length,
                      minfasta=minfasta,
                      args=args,
                      )

def binning(logger, data, minfasta,
            binned_length, contig_dict, model_path,
            output, device, environment, *, args):
    """
    Clustering the contigs to get the final bins.

    contig_dict: {contig_id:seq,...}
    recluster: if reclustering
    model_path: path to the trained model
    """
    from .cluster import cluster
    logger.info('Start binning.')

    is_combined, n_sample, data, model = binning_preprocess(data, getattr(args, 'depth_metabat2', None), model_path, environment, device)

    cluster(
        logger,
        model=model,
        data=data,
        device=device,
        is_combined=is_combined,
        n_sample=n_sample,
        out=output,
        contig_dict=contig_dict,
        args=args,
        binned_length=binned_length,
        minfasta=minfasta)


def single_easy_binning(args, logger, binned_length,
                        must_link_threshold,
                        contig_dict, output, device):
    """
    contain `generate_cannot_links`, `generate_sequence_features_single`, `train`, `bin` in one command for single-sample and co-assembly binning
    """
    logger.info('Generating training data...')
    if args.depth_metabat2 is None and args.bams is None:
        sys.stderr.write(
            f"Error: You need to input bam files if you want to calculate coverage features.\n")
        sys.exit(1)

    if args.bams is not None and args.depth_metabat2 is not None:
        logger.info('We will use abundance information from Metabat2.')

    if args.depth_metabat2 and args.environment is None:
        sys.stderr.write(
            f"Error: You need to use our provided model if you provide depth file from Metabat2.\n")
        sys.exit(1)

    generate_sequence_features_single(
        logger,
        args.contig_fasta,
        args.bams,
        binned_length,
        must_link_threshold,
        args.num_process,
        output,
        only_kmer=args.depth_metabat2)

    data_path = os.path.join(output, 'data.csv')
    if not args.depth_metabat2:
        data_split_path = os.path.join(output, 'data_split.csv')

    if args.environment is None:
        if args.training_type == 'semi':
            logger.info('Running mmseqs and generate cannot-link file.')
            predict_taxonomy(
                logger,
                args.contig_fasta,
                args.cannot_name,
                None if args.taxonomy_results_fname is None else args.taxonomy_results_fname[0],
                args.GTDB_reference,
                args.num_process,
                binned_length,
                must_link_threshold,
                output)
            logger.info('Training model and clustering.')
            training(logger, [args.contig_fasta],
                     args.num_process, [data_path], [data_split_path],
                     [os.path.join(output, 'cannot', 'cannot.txt')],
                     args.batchsize, args.epoches, output, device,
                     args.ratio, args.min_len,  mode='single',
                     orf_finder=args.orf_finder, prodigal_output_faa=args.prodigal_output_faa,
                     training_mode='semi')
        else:
            training(logger, None,
                     args.num_process, [data_path], [data_split_path],
                     None, args.batchsize, args.epoches, output, device, None, None,
                     mode='single', orf_finder=None, prodigal_output_faa=args.prodigal_output_faa, training_mode='self')

    binning_kwargs = {
        'logger': logger,
        'data': data_path,
        'args': args,
        'minfasta': args.minfasta_kb * 1000,
        'binned_length': binned_length,
        'contig_dict': contig_dict,
        'model_path':
                os.path.join(output, 'model.h5') \
                        if args.environment is None \
                        else None,
        'output': output,
        'device': device,
        'environment': args.environment,
    }

    if args.sequencing_type == 'short_read':
        binning(**binning_kwargs)
    else:
        binning_long(**binning_kwargs)


def multi_easy_binning(args, logger, output, device):
    """
    contain `generate_cannot_links`, `generate_sequence_features_multi`, `train`, `bin` in one command for multi-sample binning
    """
    logger.info('Performing multi-sample binning')
    logger.info('Generating training data...')

    sample_list = generate_sequence_features_multi(
        logger,
        args.contig_fasta,
        args.bams,
        args.num_process,
        args.separator,
        args.ratio,
        args.min_len,
        args.ml_threshold,
        output,)

    for sample_index, sample in enumerate(sample_list):
        logger.info(
            f'Running generating cannot-link file for {sample}')
        sample_fasta = os.path.join(
            output, 'samples', '{}.fa'.format(sample))
        sample_data = os.path.join(output, 'samples', sample, 'data.csv')
        sample_data_split = os.path.join(
            output, 'samples', sample, 'data_split.csv')

        binned_short, must_link_threshold, contig_dict = process_fasta(sample_fasta, args.ratio)

        if args.min_len is None:
            binned_length = 1000 if binned_short else 2500
        else:
            binned_length = args.min_len
        if args.ml_threshold is not None:
            must_link_threshold = args.ml_threshold
        logger.info('Training model and clustering for {}.'.format(sample))
        if args.training_type == 'semi':
            logger.debug(f'Running taxonomic prediction (semi-supervised mode) for {sample}')
            predict_taxonomy(
                logger,
                sample_fasta,
                sample,
                None if args.taxonomy_results_fname is None else args.taxonomy_results_fname[sample_index],
                args.GTDB_reference,
                args.num_process,
                binned_length,
                must_link_threshold,
                os.path.join(output, 'samples', sample))

            sample_cannot = os.path.join(
                output, 'samples', sample, 'cannot/{}.txt'.format(sample))
            training(logger, [sample_fasta], args.num_process,
                     [sample_data], [sample_data_split], [sample_cannot],
                     args.batchsize, args.epoches, os.path.join(output, 'samples', sample),
                     device, args.ratio, args.min_len, mode='single', orf_finder=args.orf_finder, prodigal_output_faa=args.prodigal_output_faa, training_mode='semi')
        else:
            training(logger, None, args.num_process,
                     [sample_data], [sample_data_split], None,
                     args.batchsize, args.epoches, os.path.join(output, 'samples', sample),
                     device, None, None, mode='single', orf_finder=None, prodigal_output_faa=args.prodigal_output_faa, training_mode='self')

        binning_kwargs = {
            'logger': logger,
            'args': args,
            'data': sample_data,
            'minfasta': args.minfasta_kb * 1000,
            'binned_length': binned_length,
            'contig_dict': contig_dict,
            'model_path': os.path.join(output, 'samples', sample, 'model.h5'),
            'output': os.path.join(output, 'samples', sample),
            'device': device,
            'environment': None,
        }

        if args.sequencing_type == 'short_read':
            binning(**binning_kwargs)
        else:
            binning_long(**binning_kwargs)

    os.makedirs(os.path.join(output, 'bins'), exist_ok=True)
    for sample in sample_list:
        if args.sequencing_type != 'short_read' or (not args.is_semibin2 and not args.recluster) or (args.is_semibin2 and args.recluster):
            bin_dir_name = 'output_bins'
        else:
            bin_dir_name = 'output_recluster_bins' if args.recluster else 'output_prerecluster_bins'
        for bf in os.listdir(os.path.join(output, 'samples', sample, bin_dir_name)):
            original_path = os.path.join(output, 'samples', sample, bin_dir_name, bf)
            new_file = '{0}_{1}'.format(sample, bf)
            new_path = os.path.join(output, 'bins', new_file)
            shutil.copyfile(original_path, new_path)

def main2(args=None, is_semibin2=True):
    import tempfile

    if args is None:
        args = sys.argv[1:]
    args = parse_args(args, is_semibin2)

    logger = logging.getLogger('SemiBin')
    if args.verbose:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logger.setLevel(loglevel)
    try:
        import coloredlogs
        coloredlogs.install(level=loglevel, logger=logger)
    except ImportError:
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
        logger.addHandler(sh)

    if args.cmd not in ['download_GTDB', 'check_install']:
        os.makedirs(args.output, exist_ok=True)
        fh = logging.FileHandler(os.path.join(args.output, "SemiBinRun.log"))
        fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
        logger.addHandler(fh)

    if args.verbose and args.quiet:
        logger.warning('Both verbose and quiet are set, output will be verbose')

    if sys.version_info.major < 3 or sys.version_info.minor < 7:
        logger.warning(f'You are using Python {sys.version_info.major}.{sys.version_info.minor} ({sys.version}), but SemiBin requires Python 3.7 or higher. Please upgrade your Python version.')
        logger.warning(f'If you are using conda, you can run `conda install python=3.7` to upgrade your Python version.')
        logger.warning(f'SemiBin will keep going, but it may not work properly.')

    validate_normalize_args(logger, args)
    if args.cmd == 'check_install':
        check_install(True, allow_missing_mmseqs2=args.allow_missing_mmseqs2)

    if is_semibin2 and getattr(args, 'training_type', None) == 'semi':
        logger.info('Currently using semi-supervised mode. This is generally only useful for backwards compability.')

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train', 'train_semi', 'bin', 'train_self', 'bin_long']:
        import torch
        if args.engine == 'cpu':
            device = torch.device("cpu")
            logger.info('Running with CPU.')

        else:
            if torch.cuda.is_available():
                device = torch.device("cuda")
                logger.info('Running with GPU.')
            else:
                device = torch.device("cpu")
                logger.info('Did not detect GPU, using CPU.')

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'generate_cannot_links', 'train_semi', 'train', 'bin',
                    'generate_sequence_features_single', 'generate_sequence_features_multi', 'train_self', 'bin_long']:
        tmp_output = args.tmp_output
        if tmp_output is not None:
            os.environ['TMPDIR'] = tmp_output
            os.makedirs(tmp_output, exist_ok=True)

    with tempfile.TemporaryDirectory() as tdir:
        if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'generate_sequence_features_single', 'generate_sequence_features_multi']:
            bams_new = []
            if args.bams is not None:
                for bam in args.bams:
                    bam_split = bam.split('.')
                    if bam_split[-1] == 'cram':
                        output_bam = bam.split('/')[-1][:-4] + 'bam'
                        output_bam = os.path.join(tdir, output_bam)
                        with open(output_bam, 'wb') as cram_out:
                            subprocess.check_call(
                                ['samtools', 'view',
                                 '-bS',
                                 '-@', str(args.num_process),
                                 '-T', args.contig_fasta,
                                 bam],
                                stdout=cram_out)
                        bams_new.append(output_bam)
                    else:
                        bams_new.append(bam)
                args.bams = bams_new

        if args.cmd in ['generate_cannot_links', 'generate_sequence_features_single', 'bin','single_easy_bin', 'bin_long']:
            binned_short, must_link_threshold, contig_dict = process_fasta(args.contig_fasta, args.ratio)
            if args.min_len is None:
                binned_length = 1000 if binned_short else 2500
            else:
                binned_length = args.min_len
            if not contig_dict:
                logger.error(f'Input file {args.contig_fasta} is empty. Please check inputs.')
                sys.exit(1)
            n_pass = sum(len(c) >= binned_length for c in contig_dict.values())
            if n_pass == 0:
                logger.error(f'Input file {args.contig_fasta} contains {len(contig_dict)} contigs, but all are shorter than {binned_length} basepairs.')
                sys.exit(1)
            elif n_pass < 4:
                logger.error(f'There are {len(contig_dict)} contigs in input file {args.contig_fasta}, but only {n_pass} contain(s) at least {binned_length} basepairs.')
                sys.exit(1)

        if args.cmd in ['generate_cannot_links', 'generate_sequence_features_single', 'generate_sequence_features_multi', 'single_easy_bin', 'multi_easy_bin']:
            if args.ml_threshold is not None:
                must_link_threshold = args.ml_threshold

        if args.cmd == 'download_GTDB':
            from .gtdb import find_or_download_gtdb
            find_or_download_gtdb(logger, args.GTDB_reference, args.force)

        if args.cmd == 'generate_cannot_links':
            predict_taxonomy(
                logger,
                args.contig_fasta,
                args.cannot_name,
                taxonomy_results_fname=None if args.taxonomy_results_fname is None else args.taxonomy_results_fname[0],
                GTDB_reference=args.GTDB_reference,
                num_process=args.num_process,
                binned_length=binned_length,
                must_link_threshold=must_link_threshold,
                output=args.output)

        if args.cmd == 'generate_sequence_features_single':
            generate_sequence_features_single(
                logger,
                args.contig_fasta,
                args.bams,
                binned_length,
                must_link_threshold,
                args.num_process,
                args.output,
                args.kmer)

        if args.cmd == 'generate_sequence_features_multi':
            generate_sequence_features_multi(
                logger,
                args.contig_fasta,
                args.bams,
                args.num_process,
                args.separator,
                args.ratio,
                args.min_len,
                args.ml_threshold,
                args.output)

        if args.cmd in ['train', 'train_semi']:
            if args.random_seed is not None:
                set_random_seed(args.random_seed)
            training(logger, args.contig_fasta, args.num_process,
                     args.data, args.data_split, args.cannot_link,
                     args.batchsize, args.epoches, args.output, device, args.ratio, args.min_len, args.mode, orf_finder=args.orf_finder, training_mode='semi')

        if args.cmd == 'train_self':
            if args.random_seed is not None:
                set_random_seed(args.random_seed)
            training(logger, None, args.num_process,
                     args.data, args.data_split, None,
                     args.batchsize, args.epoches, args.output, device, None, None,
                     mode=args.mode, orf_finder=None,
                     prodigal_output_faa=None, training_mode='self')


        if args.cmd == 'bin':
            if args.random_seed is not None:
                set_random_seed(args.random_seed)
            binning(logger, args.data, args.minfasta_kb * 1000, binned_length,
                    environment=args.environment, contig_dict=contig_dict,
                    model_path=args.model_path, output=args.output,
                    device=device, args=args)

        if args.cmd == 'bin_long':
            if args.random_seed is not None:
                set_random_seed(args.random_seed)
            binning_long(logger, args.data, args.minfasta_kb * 1000,
                    binned_length, environment=args.environment,
                    contig_dict=contig_dict, model_path=args.model_path,
                    output=args.output, device=device, args=args)

        if args.cmd == 'single_easy_bin':
            check_install(False, args.orf_finder, allow_missing_mmseqs2=(args.environment is not None or args.training_type == 'self'))
            if args.random_seed is not None:
                set_random_seed(args.random_seed)
            if args.environment is not None:
                if args.depth_metabat2 is None:
                    if len(args.bams) != 1:
                        sys.stderr.write(
                            f"Error: provided pretrained model only used in single-sample binning!\n")
                        sys.exit(1)

            single_easy_binning(
                args,
                logger,
                binned_length,
                must_link_threshold,
                contig_dict,
                args.output,
                device)

        if args.cmd == 'multi_easy_bin':
            check_install(False, args.orf_finder, args.training_type == 'self')
            if args.random_seed is not None:
                set_random_seed(args.random_seed)
            multi_easy_binning(
                args,
                logger,
                args.output,
                device)

        if args.cmd == 'concatenate_fasta':
            from .utils import concatenate_fasta
            ofname = concatenate_fasta(args.contig_fasta, args.min_length, args.output, args.separator, args.output_compression)
            logger.info(f'Concatenated contigs written to {ofname}')

        print('''If you find SemiBin useful, please cite:
        Pan, S.; Zhu, C.; Zhao, XM.; Coelho, LP. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. Nat Commun 13, 2326 (2022). https://doi.org/10.1038/s41467-022-29843-y

        Pan, S.; Zhao, XM; Coelho, LP. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. Bioinformatics Volume 39, Issue Supplement_1, June 2023, Pages i21–i29. https://doi.org/10.1093/bioinformatics/btad209

''')


def main1(args=None):
    main2(args, is_semibin2=False)

if __name__ == '__main__':
    main2()
