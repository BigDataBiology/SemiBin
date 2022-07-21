from torch.serialization import SourceChangeWarning
import warnings
warnings.filterwarnings("ignore", category=SourceChangeWarning)
from .semibin_version import __version__ as ver
import argparse
import logging
import os
import subprocess
from atomicwrites import atomic_write
import shutil
import sys
from itertools import groupby

from .utils import validate_normalize_args, get_must_link_threshold, generate_cannot_link, \
    set_random_seed, process_fasta, split_data, get_model_path
from .gtdb import find_or_download_gtdb
from .generate_coverage import generate_cov, combine_cov
from .generate_kmer import generate_kmer_features_from_fasta
from .cluster import cluster
from .fasta import fasta_iter
from .error import LoggingPool


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural '
                                                 'network for metagenomic binning')

    parser.version = ver

    parser.add_argument('-v',
                        '--version',
                        action='version',
                        help='Print the version number')

    subparsers = parser.add_subparsers(title='SemiBin subcommands',
                                       dest='cmd',
                                       metavar='')

    single_easy_bin = subparsers.add_parser('single_easy_bin',
                                            help='Bin contigs (single or co-assembly) using one command.')

    multi_easy_bin = subparsers.add_parser('multi_easy_bin',
                                            help='Bin contigs (multi-sample mode) using one command.')

    generate_cannot_links = subparsers.add_parser('generate_cannot_links', aliases=['predict_taxonomy',],
                                             help='Run the contig annotation using mmseqs '
                                                  'with GTDB reference genome and generate '
                                                  'cannot-link file used in the semi-supervised deep learning model training. '
                                                  'This will download the GTDB database if not downloaded before.')

    generate_sequence_features_single = subparsers.add_parser('generate_sequence_features_single', aliases=['generate_sequence_features_single'],
                                            help='Generate sequence features (kmer and abundance) as training data'
                                                  ' for semi-supervised deep learning model training (single or co-assembly mode).'
                                                  ' This will produce the data.csv and data_split.csv files.'
                                                  )


    generate_sequence_features_multi = subparsers.add_parser('generate_sequence_features_multi', aliases=['generate_sequence_features_multi'],
                                            help='Generate sequence features (kmer and abundance) as training data'
                                                  ' for semi-supervised deep learning model training (multi-sample mode).'
                                                  ' This will produce the data.csv and data_split.csv files.'
                                                  )

    download_GTDB = subparsers.add_parser('download_GTDB', help='Download GTDB reference genomes.')

    check_install = subparsers.add_parser('check_install', help = 'Check whether required dependencies are present.')

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

    concatenate_fasta.add_argument('-m',
                        required=False,
                        type=int,
                        help='Discard sequences below this length (default:0)',
                        default=0,
                        dest='min_length')

    training = subparsers.add_parser('train',
                                    help='Train the model.')

    binning = subparsers.add_parser('bin',
                                    help='Group the contigs into bins.')

    download_GTDB.add_argument('-f', '--force',
                            required=False,
                            help='Redownload GTDB even if files are found',
                            dest='force',
                            action='store_true',
                            default=None)

    training.add_argument('--data',
                         required=True,
                         nargs='*',
                         help='Path to the input data.csv file.',
                         dest='data',
                         default=None,
                         )
    training.add_argument('--data-split',
                         required=True,
                         nargs='*',
                         help='Path to the input data_split.csv file.',
                         dest='data_split',
                         default=None,
                         )
    training.add_argument('-c', '--cannot-link',
                         required=True,
                          nargs='*',
                         help='Path to the input cannot link file. '
                         'The file format: `contig_1,contig_2` '
                         '(one row for each cannot link constraint).',
                         dest='cannot_link',
                         default=None)

    training.add_argument('--epoches',
                   required=False,
                   type=int,
                   help='Number of epoches used in the training process (Default: 20).',
                   dest='epoches',
                   default=20)

    training.add_argument('--batch-size',
                   required=False,
                   type=int,
                   help='Batch size used in the training process (Default: 2048).',
                   dest='batchsize',
                   default=2048, )

    training.add_argument('--mode',
                          required=True,
                          type=str,
                          help='[single/several] Train models from one (single) or more samples (several). '
                          'In `several` mode, you must provide data, data_split, cannot, and fasta files for corresponding samples in the same order. '
                          'Note: You can only use `several` mode when performing single-sample binning. Training from several samples with multi-sample binning is not supported.',
                          dest='mode',
                          default='single')

    training.add_argument('-i', '--input-fasta',
                   required=True,
                   nargs='*',
                   help='Path to the input fasta file.',
                   dest='contig_fasta',
                   default=None, )

    training.add_argument('-o', '--output',
                   required=True,
                   help='Output directory (will be created if non-existent)',
                   dest='output',
                   default=None,
                   )

    binning.add_argument('--data',
                         required=True,
                         help='Path to the input data.csv file.',
                         dest='data',
                         default=None,)

    binning.add_argument('--max-edges',
                   required=False,
                   type=int,
                   help='The maximum number of edges that can be connected to one contig (Default: 200).',
                   dest='max_edges',
                   default=200)
    binning.add_argument('--max-node',
                   required=False,
                   type=float,
                   dest='max_node',
                   default=1,
                   help='Fraction of contigs that considered to be binned (should be between 0 and 1; default: 1).')

    binning.add_argument('--model',
                         required=False,
                         type=str,
                         dest='model_path',
                         default=None,
                         help='Path to the trained semi-supervised deep learning model.')

    for p in [single_easy_bin, multi_easy_bin, training, binning]:
        p.add_argument('--orf-finder',
                       required=False,
                       type=str,
                       help='ORF finder used to estimate the number of bins (prodigal/fraggenescan)',
                       dest='orf_finder',
                       default='prodigal')

    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links, training, binning]:
        p.add_argument('--tmpdir',
                       required=False,
                       type=str,
                       help='option to set temporary directory',
                       dest='tmp_output',
                       default=None,
    )

    for p in [training, generate_cannot_links, binning, single_easy_bin, multi_easy_bin, generate_sequence_features_single, generate_sequence_features_multi]:
        p.add_argument('-p', '--processes', '-t', '--threads',
                   required=False,
                   type=int,
                   help='Number of CPUs used (pass the value 0 to use all CPUs, default: 0)',
                   dest='num_process',
                   default=0,
                   metavar=''
                   )

    for p in [single_easy_bin, binning]:
        p.add_argument('--environment',
                       required=False,
                       help='Environment for the built-in model (available choices: human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global).',
                       dest='environment',
                       default=None,
                       )

    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links, generate_sequence_features_single, generate_sequence_features_multi, binning]:
        p.add_argument('-i', '--input-fasta',
                                required=True,
                                help='Path to the input fasta file.',
                                dest='contig_fasta',
                                default=None,)
        p.add_argument('-o', '--output',
                            required=True,
                            help='Output directory (will be created if non-existent)',
                            dest='output',
                            default=None,
                            )
    for p in [single_easy_bin, multi_easy_bin, generate_cannot_links, generate_sequence_features_single, generate_sequence_features_multi, binning, training]:
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

    for p in [single_easy_bin, multi_easy_bin, generate_sequence_features_single, generate_sequence_features_multi]:
        p.add_argument('-b', '--input-bam',
                              required=True,
                              nargs='*',
                              help='Path to the input BAM file. '
                                   'If using multiple samples, you can input multiple files.',
                              dest='bams',
                              default=None,
                              )

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
                            help='Name for the cannot-link file(default: cannot).',
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

    for p in [binning, single_easy_bin, multi_easy_bin]:
        p.add_argument('--minfasta-kbs',
                            required=False,
                            type=int,
                            help='minimum bin size in Kbps (Default: 200).',
                            dest='minfasta_kb',
                            default=200,
                            metavar='')

        p.add_argument('--no-recluster',
                           required=False,
                           help='Do not recluster bins.',
                           dest='no_recluster',
                           action='store_true', )
        p.add_argument('--recluster',
                           required=False,
                           help='[Deprecated] Does nothing (current default is to perform clustering)',
                           dest='recluster',
                           action='store_true', )
    for p in [single_easy_bin, multi_easy_bin]:

        p.add_argument('--epoches',
                          required=False,
                          type=int,
                          help='Number of epoches used in the training process (Default: 20).',
                          dest='epoches',
                          default=20)

        p.add_argument('--batch-size',
                          required=False,
                          type=int,
                          help='Batch size used in the training process (Default: 2048).',
                          dest='batchsize',
                          default=2048,)

        p.add_argument('--max-edges',
                          required=False,
                          type=int,
                          help='The maximum number of edges that can be connected to one contig (Default: 200).',
                          dest='max_edges',
                          default=200)

        p.add_argument('--max-node',
                          required=False,
                          type=float,
                          dest='max_node',
                          default=1,
                          help='Fraction of contigs that considered to be binned (should be between 0 and 1; default: 1).')

    for p in [multi_easy_bin, generate_sequence_features_multi, concatenate_fasta]:
        p.add_argument('-s', '--separator',
                           required=False,
                           type=str,
                           help='Used when multiple samples binning to separate sample name and contig name (Default is :).',
                           dest='separator',
                           default=':',
                           metavar='')

    for p in [training, binning, single_easy_bin, multi_easy_bin]:
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

    for p in [single_easy_bin, multi_easy_bin, training, binning]:
        p.add_argument('--engine',
                       required=False,
                       type=str,
                       help='device used to train the model (auto/gpu/cpu, auto means if SemiBin detects the gpu, SemiBin will use GPU)',
                       dest='engine',
                       default='auto')


    if not args:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args(args)
    if hasattr(args, 'no_recluster'):
        args.recluster = not args.no_recluster
    return args


def _checkback(msg):
    msg[1].info('Processed:{}'.format(msg[0]))

def check_install(verbose, orf_finder=None):
    '''Executes check_install subcommand, which checks for dependencies

    Parameters
    ----------
    verbose : boolean. If true, then prints out all the paths
    orf_finder : str, optional
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
            if dep not in ['FragGeneScan', 'prodigal']:
                sys.stderr.write(
                    f"Error: {dep} does not seem to be installed!\n")
                missing_deps = True
            elif dep == 'prodigal':
                if not has_fgs:
                    sys.stderr.write(
                            'Error: neither prodigal nor FragGeneScan appear to be available!\n'
                            'At least one of them is necessary to run SemiBin\n')
                    missing_deps = True
                else:
                    if verbose or orf_finder != 'fraggenescan':
                        sys.stderr.write(
                            'Warning: prodigal does not appear to be available. You must use the `--orf-finder fraggenescan` option.\n')
                    missing_deps = True
        else:
            if dep == 'FragGeneScan':
                has_fgs = True
            if verbose:
                print(f'\t{dep:16}: {p}')
    if missing_deps:
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
                         must_link_threshold, num_process, output):
    """
    Generate data.csv and data_split.csv for training and clustering of single-sample and co-assembly binning mode.
    data.csv has the features(kmer and abundance) for original contigs.
    data_split.csv has the features(kmer and abundace) for contigs that are breaked up as must-link pair.

    """
    import pandas as pd
    n_sample = len(bams)
    is_combined = n_sample >= 5
    bam_list = bams
    pool = LoggingPool(num_process) if num_process != 0 else LoggingPool()

    logger.info('Calculating coverage for every sample.')

    for bam_index in range(n_sample):
        pool.apply_async(
            generate_cov,
            args=(
                bam_list[bam_index],
                bam_index,
                output,
                must_link_threshold,
                is_combined,
                binned_length,
                logger,
                None
            ),
            callback=_checkback)
    pool.close()
    pool.join()

    for bam_index, bam_file in enumerate(bam_list):
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

    logger.info('Start generating kmer features from fasta file.')
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
    if num_process != 0:
        pool = LoggingPool(num_process)
    else:
        pool = LoggingPool()

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

    for bam_index in range(n_sample):
        pool.apply_async(generate_cov,
                         args=(
                             bam_list[bam_index],
                             bam_index,
                             os.path.join(output, 'samples'),
                             must_link_threshold,
                             is_combined,
                             binning_threshold,
                             logger,
                             separator,
                         ),
                         callback=_checkback)
    pool.close()
    pool.join()

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
             epoches,  output, device, ratio, min_length, mode, orf_finder = 'prodigal'):
    """
    Training the model

    model: [single/several]
    """
    from .semi_supervised_model import train
    import pandas as pd
    binned_lengths = []

    if mode == 'single':
        logger.info('Start training from one sample.')
        if min_length is None:
            binned_short, _, _ = process_fasta(contig_fasta[0], ratio)
            binned_lengths.append(1000) if binned_short else binned_lengths.append(2500)
        else:
            binned_lengths.append(min_length)
        data_ = pd.read_csv(data[0], index_col=0)
        col_name = data_.columns.tolist()[-1].split('_')[-1]
        is_combined = False if col_name == 'var' else True

    else:
        logger.info('Start training from multiple samples.')
        is_combined = False
        for contig_index in contig_fasta:
            if min_length is None:
                binned_short, _, _ = process_fasta(contig_index, ratio)
                binned_lengths.append(1000) if binned_short else binned_lengths.append(2500)
            else:
                binned_lengths.append(min_length)

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
        orf_finder=orf_finder)


def binning(logger, num_process, data,
            max_edges, max_node, minfasta,
            binned_length, contig_dict, recluster,model_path,
            random_seed,output, device, environment, orf_finder = 'prodigal'):
    """
    Clustering the contigs to get the final bins.

    contig_dict: {contig_id:seq,...}
    recluster: if reclustering
    model_path: path to the trained model
    """
    import torch
    import pandas as pd
    logger.info('Start binning.')

    data = pd.read_csv(data, index_col=0)
    data.index = data.index.astype(str)

    col_name = data.columns.tolist()[-1].split('_')[-1]
    is_combined = False if col_name == 'var' else True
    n_sample = (data.shape[1] - 136) // 2 if not is_combined else (data.shape[1] - 136)
    model_path = model_path if environment is None else get_model_path(environment)
    if environment is not None:
        if data.shape[1] != 138:
            sys.stderr.write(f"Error: provided pretrained model only used in single-sample binning!\n")
            sys.exit(1)

    if device == torch.device('cpu'):
        model = torch.load(model_path, map_location=torch.device('cpu'))
    else:
        model = torch.load(model_path)
    cluster(
        model,
        data,
        device,
        max_edges,
        max_node,
        is_combined,
        logger,
        n_sample,
        output,
        contig_dict,
        binned_length,
        num_process,
        minfasta,
        recluster,
        random_seed,
        orf_finder=orf_finder)


def single_easy_binning(args, logger, binned_length,
                        must_link_threshold,
                        contig_dict, recluster,random_seed, output, device, environment, orf_finder = 'prodigal'):
    """
    contain `generate_cannot_links`, `generate_sequence_features_single`, `train`, `bin` in one command for single-sample and co-assembly binning
    """
    logger.info('Generate training data.')
    generate_sequence_features_single(
        logger,
        args.contig_fasta,
        args.bams,
        binned_length,
        must_link_threshold,
        args.num_process,
        output)
    data_path = os.path.join(output, 'data.csv')
    data_split_path = os.path.join(output, 'data_split.csv')

    if environment is None:
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
                 args.batchsize, args.epoches, output, device, args.ratio, args.min_len,  mode='single', orf_finder=orf_finder)

    binning(logger, args.num_process, data_path,
            args.max_edges, args.max_node, args.minfasta_kb * 1000,
            binned_length, contig_dict, recluster,
            os.path.join(output, 'model.h5') if environment is None else None,
            random_seed, output,  device, environment if environment is not None else None, orf_finder=orf_finder)


def multi_easy_binning(args, logger, recluster,
                       random_seed, output, device, orf_finder='prodigal'):
    """
    contain `generate_cannot_links`, `generate_sequence_features_multi`, `train`, `bin` in one command for multi-sample binning
    """
    logger.info('Multi-sample binning.')
    logger.info('Generate training data.')

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
            'Running mmseqs and generate cannot-link file of {}.'.format(sample))
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
        logger.info('Training model and clustering for {}.'.format(sample))
        training(logger, [sample_fasta], args.num_process,
                 [sample_data], [sample_data_split], [sample_cannot],
                 args.batchsize, args.epoches, os.path.join(output, 'samples', sample),
                 device, args.ratio, args.min_len, mode='single', orf_finder=orf_finder)

        binning(logger, args.num_process, sample_data,
                args.max_edges, args.max_node, args.minfasta_kb * 1000,
                binned_length, contig_dict, recluster,
                os.path.join(output, 'samples', sample, 'model.h5'), random_seed ,
                os.path.join(output, 'samples', sample),  device, None, orf_finder=orf_finder)

    os.makedirs(os.path.join(output, 'bins'), exist_ok=True)
    for sample in sample_list:
        bin_dir_name = 'output_recluster_bins' if recluster else 'output_bins'
        for bf in os.listdir(os.path.join(output, 'samples', sample, bin_dir_name)):
            original_path = os.path.join(output, 'samples', sample, bin_dir_name, bf)
            new_file = '{0}_{1}'.format(sample, bf)
            new_path = os.path.join(output, 'bins', new_file)
            shutil.copyfile(original_path, new_path)


def main():
    args = sys.argv[1:]
    args = parse_args(args)

    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    validate_normalize_args(logger, args)

    if args.cmd == 'check_install':
        check_install(True)

    if args.cmd not in ['download_GTDB', 'check_install']:
        out = args.output
        os.makedirs(out, exist_ok=True)

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train', 'bin']:
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
                logger.info('Do not detect GPU. Running with CPU.')

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'generate_cannot_links', 'train', 'bin']:
        tmp_output = args.tmp_output
        if tmp_output is not None:
            os.environ['TMPDIR'] = tmp_output
            os.makedirs(tmp_output, exist_ok=True)

    if args.cmd in ['generate_cannot_links', 'generate_sequence_features_single', 'bin','single_easy_bin']:
        binned_short, must_link_threshold, contig_dict = process_fasta(args.contig_fasta, args.ratio)
        if args.min_len is None:
            binned_length = 1000 if binned_short else 2500
        else:
            binned_length = args.min_len

    if args.cmd in ['generate_cannot_links', 'generate_sequence_features_single', 'generate_sequence_features_multi', 'single_easy_bin', 'multi_easy_bin']:
        if args.ml_threshold is not None:
            must_link_threshold = args.ml_threshold

    if args.cmd == 'download_GTDB':
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
            output=out)

    if args.cmd == 'generate_sequence_features_single':
        generate_sequence_features_single(
            logger,
            args.contig_fasta,
            args.bams,
            binned_length,
            must_link_threshold,
            args.num_process,
            out)

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
            out)

    if args.cmd == 'train':
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        training(logger, args.contig_fasta, args.num_process,
                 args.data, args.data_split, args.cannot_link,
                 args.batchsize, args.epoches, out, device, args.ratio, args.min_len, args.mode, orf_finder=args.orf_finder)


    if args.cmd == 'bin':
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        binning(logger, args.num_process, args.data, args.max_edges,
                args.max_node, args.minfasta_kb * 1000, binned_length,
                contig_dict, args.recluster, args.model_path, args.random_seed,out, device, args.environment, orf_finder=args.orf_finder)


    if args.cmd == 'single_easy_bin':
        check_install(False, args.orf_finder)
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        if args.environment is not None:
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
            args.recluster,
            args.random_seed,
            out,
            device,
            args.environment,
            orf_finder=args.orf_finder)

    if args.cmd == 'multi_easy_bin':
        check_install(False, args.orf_finder)
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        multi_easy_binning(
            args,
            logger,
            args.recluster,
            args.random_seed,
            out,
            device,
            orf_finder=args.orf_finder)

    if args.cmd == 'concatenate_fasta':
        from .utils import concatenate_fasta
        concatenate_fasta(args.contig_fasta, args.min_length, out, args.separator)

    print('If you find SemiBin useful, please cite:\n Pan, S., Zhu, C., Zhao, XM. et al. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. Nat Commun 13, 2326 (2022). https://doi.org/10.1038/s41467-022-29843-y.')

if __name__ == '__main__':
    main()
