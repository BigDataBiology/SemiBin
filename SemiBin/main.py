from torch.serialization import SourceChangeWarning
import warnings
warnings.filterwarnings("ignore", category=SourceChangeWarning)
from .semibin_version import __version__ as ver
import argparse
import logging
import os
import multiprocessing
import subprocess
from atomicwrites import atomic_write
import shutil
import sys
from Bio import SeqIO
from .utils import validate_args, get_threshold, generate_cannot_link, \
    download, set_random_seed, unzip_fasta, process_fasta, split_data, get_model_path
from .generate_coverage import generate_cov, combine_cov
from .generate_kmer import generate_kmer_features_from_fasta
from Bio.SeqRecord import SeqRecord
from .semi_supervised_model import train
from .cluster import cluster
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

    predict_taxonomy = subparsers.add_parser('predict_taxonomy',
                                             help='Run the contig annotation using mmseqs '
                                                  'with GTDB reference genome and generate '
                                                  'cannot-link file used in the semi-supervised deep learning model training. '
                                                  'This will download the GTDB database if not downloaded before.')

    generate_data_single = subparsers.add_parser('generate_data_single',
                                            help='Generate training data (files data.csv and data_split.csv) '
                                                  'for semi-supervised deep learning model training (single or co-assembly).')


    generate_data_multi = subparsers.add_parser('generate_data_multi',
                                            help='Generate training data(files data.csv and data_split.csv) '
                                                'for the semi-supervised deep learning model training (multi-sample)')

    download_GTDB = subparsers.add_parser('download_GTDB', help='Download GTDB reference genomes.')

    training = subparsers.add_parser('train',
                                    help='Train the model.')

    binning = subparsers.add_parser('bin',
                                    help='Binning the contigs into bins.')

    download_GTDB.add_argument('-r', '--reference-db',
                            required=False,
                            help='GTDB reference file path to download(~/path). (Default: $HOME/.cache/SemiBin/mmseqs2-GTDB). If not set --reference-db, we will download GTDB to the default path.',
                            dest='GTDB_reference',
                            metavar='',
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
                          help='[single/several] Train models from one sample or several samples (train model across several samples with single-sample binning can get better pre-trained model). '
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

    binning.add_argument('--minfasta-kbs',
                            required=False,
                            type=int,
                            help='minimum bin size in Kbps (Default: 200).',
                            dest='minfasta_kb',
                            default=200,
                            metavar='')
    binning.add_argument('--recluster',
                   required=False,
                   help='recluster bins.',
                   dest='recluster',
                   action='store_true', )
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

    training.add_argument('-p', '--processes', '-t', '--threads',
                   required=False,
                   type=int,
                   help='Number of CPUs used (pass the value 0 to use all CPUs, default: 0)',
                   dest='num_process',
                   default=0,
                   metavar=''
                   )
    training.add_argument('-b', '--input-bam',
                   required=False,
                   nargs='*',
                   help='Path to the input BAM file(Only need when training mode is single). '
                        'If using multiple samples, you can input multiple files.',
                   dest='bams',
                   default=None,
                   )
    predict_taxonomy.add_argument('-p', '--processes', '-t', '--threads',
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
                       help='environment for the built-in model(human_gut/dog_gut/ocean).',
                       dest='environment',
                       default=None,
                       )

    for p in [single_easy_bin, multi_easy_bin, predict_taxonomy, generate_data_single, generate_data_multi, binning]:
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
    for p in [single_easy_bin, multi_easy_bin, predict_taxonomy, generate_data_single, generate_data_multi, binning, training]:
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

    for p in [single_easy_bin, multi_easy_bin, generate_data_single, generate_data_multi, binning]:
        p.add_argument('-p', '--processes', '-t', '--threads',
                                     required=False,
                                     type=int,
                                     help='Number of CPUs used (pass the value 0 to use all CPUs, default: 0)',
                                     dest='num_process',
                                     default=0,
                                     metavar=''
                                     )
        p.add_argument('-b', '--input-bam',
                              required=True,
                              nargs='*',
                              help='Path to the input BAM file. '
                                   'If using multiple samples, you can input multiple files.',
                              dest='bams',
                              default=None,
                              )

    for p in [single_easy_bin, multi_easy_bin, predict_taxonomy]:
        p.add_argument('-r', '--reference-db',
                            required=False,
                            help='GTDB reference storage path. (Default: $HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB).'
                            'If not set --reference-db and SemiBin cannot find GTDB in $HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB, SemiBin will download GTDB (Note that >100GB of disk space are required).',
                            dest='GTDB_reference',
                            metavar='',
                            default=None)

    for p in [single_easy_bin, predict_taxonomy]:
        p.add_argument('--cannot-name',
                            required=False,
                            help='Name for the cannot-link file(default: cannot).',
                            dest='cannot_name',
                            default='cannot',
                            metavar=''
                            )

    for p in [single_easy_bin, multi_easy_bin]:
        p.add_argument('--minfasta-kbs',
                            required=False,
                            type=int,
                            help='minimum bin size in Kbps (Default: 200).',
                            dest='minfasta_kb',
                            default=200,
                            metavar='')

        p.add_argument('--recluster',
                            required=False,
                            help='recluster bins.',
                            dest='recluster',
                            action ='store_true',)

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

    for p in [multi_easy_bin, generate_data_multi]:
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
                       help='Random seed to reproduce results.',
                       dest='random_seed',
                       default=None,
                       )

    for p in [predict_taxonomy, generate_data_single, generate_data_multi, single_easy_bin, multi_easy_bin]:
        p.add_argument('--ml-threshold',
                       required=False,
                       type=int,
                       help='Length threshold for generating must-link constraints. (By default, the threshold is calculated from the contig, and the default minimum value is 4,000 bp)',
                       dest='ml_threshold',
                       default=None)

    if not args:
        parser.print_help(sys.stderr)
        sys.exit()
    return parser.parse_args(args)


def _checkback(msg):
    msg[1].info('Processed:{}'.format(msg[0]))


def download_GTDB(logger,GTDB_reference):
    GTDB_default = os.path.join(
        os.environ['HOME'],
        '.cache',
        'SemiBin',
        'mmseqs2-GTDB')

    GTDB_path = GTDB_reference if GTDB_reference is not None else GTDB_default
    download(logger, GTDB_path)


def predict_taxonomy(logger, contig_fasta,
                     cannot_name, GTDB_reference, num_process,
                     binned_length, must_link_threshold,output):
    """
    Predict taxonomy using mmseqs and generate cannot-link file

    contig_fasta: contig file used
    GTDB_reference: GTDB path
    cannot_name: name for cannot-link constraint
    binned_short: threshold for contigs used in binning
    must_link_threshold: threshold of contigs for must-link constraints
    """
    import tempfile
    num_cpu = multiprocessing.cpu_count() if num_process == 0 else num_process
    GTDB_default = os.path.join(
        os.environ['HOME'],
        '.cache',
        'SemiBin',
        'mmseqs2-GTDB',
        'GTDB')
    if GTDB_reference is None:
        if not os.path.exists(GTDB_default):
            download(logger, GTDB_default)
        GTDB_reference = GTDB_default

    filtered_fasta = []
    for seq_record in SeqIO.parse(contig_fasta, "fasta"):
        if len(seq_record) > binned_length:
            filtered_fasta.append(seq_record)

    SeqIO.write(filtered_fasta, os.path.join(output, 'filtered.fasta'), 'fasta')
    filtered_fasta = os.path.join(output, 'filtered.fasta')
    subprocess.check_call(
        ['mmseqs',
         'createdb',
         filtered_fasta,
         os.path.join(output, 'contig_DB')],
        stdout=None,
    )
    if os.path.exists(os.path.join(output, 'mmseqs_contig_annotation')):
        shutil.rmtree(os.path.join(output, 'mmseqs_contig_annotation'))
    os.makedirs(os.path.join(output, 'mmseqs_contig_annotation'))
    with tempfile.TemporaryDirectory() as tdir:
        subprocess.run(
            ['mmseqs',
             'taxonomy',
             os.path.join(output, 'contig_DB'),
             GTDB_reference,
             os.path.join(output, 'mmseqs_contig_annotation/mmseqs_contig_annotation'),
             tdir,
             '--tax-lineage', '1',
             '--threads', str(num_cpu),
             ],
            check=True,
            stdout=None,
        )
    subprocess.check_call(
        ['mmseqs',
         'createtsv',
         os.path.join(output, 'contig_DB'),
         os.path.join(output, 'mmseqs_contig_annotation/mmseqs_contig_annotation'),
         os.path.join(output, 'mmseqs_contig_annotation/taxonomyResult.tsv')
         ],
        stdout=None,
    )

    namelist = []
    num_must_link = 0

    for seq_record in SeqIO.parse(filtered_fasta, "fasta"):
        if len(seq_record) > binned_length:
            namelist.append(seq_record.id)
        if len(seq_record) >= must_link_threshold:
            num_must_link += 1
    os.makedirs(os.path.join(output, 'cannot'), exist_ok=True)
    generate_cannot_link(os.path.join(output, 'mmseqs_contig_annotation/taxonomyResult.tsv'),
        namelist, num_must_link, os.path.join(output, 'cannot'), cannot_name)


def generate_data_single(logger, contig_fasta,
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

    logger.info('Start generating kmer features from fasta file.')
    kmer_whole = generate_kmer_features_from_fasta(
        contig_fasta, binned_length, 4)
    kmer_split = generate_kmer_features_from_fasta(
        contig_fasta, 1000, 4, split=True, threshold=must_link_threshold)

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


def generate_data_multi(logger, contig_fasta,
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
    from collections import defaultdict
    sample_list = list()
    contig_sample_list = []
    contig_length_list = []
    flag_name = None

    os.makedirs(os.path.join(output, 'samples'), exist_ok=True)

    for seq_record in SeqIO.parse(contig_fasta, "fasta"):
        sample_name, contig_name = seq_record.id.split(separator)
        if flag_name is None:
            flag_name = sample_name
        if sample_name == flag_name:
            rec = SeqRecord(seq_record.seq, id=contig_name, description='')
            contig_sample_list.append(rec)
        if sample_name != flag_name:
            SeqIO.write(contig_sample_list,
                        os.path.join(
                            output, 'samples', '{}.fasta'.format(flag_name)), 'fasta')
            flag_name = sample_name
            contig_sample_list = []
            rec = SeqRecord(seq_record.seq, id=contig_name, description='')
            contig_sample_list.append(rec)
        if sample_name not in sample_list:
            sample_list.append(sample_name)
        contig_length_list.append(len(seq_record))

    if contig_sample_list != []:
        SeqIO.write(contig_sample_list,
                    os.path.join(
                        output, 'samples', '{}.fasta'.format(flag_name)), 'fasta')


    must_link_threshold = get_threshold(contig_length_list) if ml_threshold is None else ml_threshold

    logger.info('Calculating coverage for every sample.')

    binning_threshold = {}
    for sample in sample_list:
        if min_length is not None:
            binning_threshold[sample] = min_length
        else:
            binned_short ,_ ,_ ,_ = process_fasta(os.path.join(output, 'samples/{}.fasta'.format(sample)), ratio)
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

        part_data = split_data(data_cov, sample, separator)
        part_data.to_csv(os.path.join(output_path, 'data_cov.csv'))

        if is_combined:
            part_data = split_data(data_split_cov, sample, separator)
            part_data.to_csv(os.path.join(
                output_path, 'data_split_cov.csv'))

        sample_contig_fasta = os.path.join(
            output, 'samples/{}.fasta'.format(sample))
        binned_length =  binning_threshold[sample]
        kmer_whole = generate_kmer_features_from_fasta(
            sample_contig_fasta, binned_length, 4)
        kmer_split = generate_kmer_features_from_fasta(
            sample_contig_fasta, 1000, 4, split=True, threshold=must_link_threshold)


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


def training(logger, contig_fasta, bams, num_process,
             data, data_split, cannot_link, batchsize,
             epoches,  output, device, ratio, min_length, mode):
    """
    Training the model

    model: [single/several]
    """

    binned_lengths = []
    num_cpu = multiprocessing.cpu_count() if num_process == 0 else num_process

    if mode == 'single':
        logger.info('Start training from one sample.')
        if min_length is None:
            binned_short, _, _, _ = process_fasta(contig_fasta[0], ratio)
            binned_lengths.append(1000) if binned_short else binned_lengths.append(2500)
        else:
            binned_lengths.append(min_length)
        n_sample = len(bams)
        is_combined = n_sample >= 5

    else:
        logger.info('Start training from multiple samples.')
        is_combined = False
        for contig_index in contig_fasta:
            if min_length is None:
                binned_short, _, _, _ = process_fasta(contig_index, ratio)
                binned_lengths.append(1000) if binned_short else binned_lengths.append(2500)
            else:
                binned_lengths.append(min_length)

    model = train(
        output,
        contig_fasta,
        binned_lengths,
        logger,
        data,
        data_split,
        cannot_link,
        is_combined,
        batchsize,
        epoches,
        device,
        num_cpu,
        mode = mode)


def binning(logger,bams, num_process, data,
            max_edges, max_node, minfasta,
            binned_length, contig_length_dict,
            contig_dict,recluster,model_path,
            random_seed,output, device, environment):
    """
    Clustering the contigs to get the final bins.

    contig_length_dict: {contig_id:length,...}
    contig_dict: {contig_id:seq,...}
    recluster: if reclustering
    model_path: path to the trained model
    """
    import torch
    import pandas as pd
    logger.info('Start binning.')
    n_sample = len(bams)
    is_combined = n_sample >= 5
    num_cpu = multiprocessing.cpu_count() if num_process == 0 else num_process
    data = pd.read_csv(data, index_col=0)
    data.index = data.index.astype(str)

    model_path = model_path if environment is None else get_model_path(environment)

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
        contig_length_dict,
        output,
        contig_dict,
        binned_length,
        num_cpu,
        minfasta,
        recluster,
        random_seed)


def single_easy_binning(args, logger, binned_length,
                        must_link_threshold, contig_length_dict,
                        contig_dict, recluster,random_seed, output, device, environment):
    """
    contain `predict_taxonomy`, `generate_data_single`, `train`, `bin` in one command for single-sample and co-assembly binning
    """
    logger.info('Generate training data.')
    generate_data_single(
        logger,
        args.contig_fasta,
        args.bams,
        binned_length,
        must_link_threshold,
        args.num_process,
        output)
    data_path = os.path.join(output, 'data.csv')
    data_split_path = os.path.join(output, 'data_split.csv')

    num_cpu = multiprocessing.cpu_count() if args.num_process == 0 else args.num_process

    if environment is None:
        logger.info('Running mmseqs and generate cannot-link file.')
        predict_taxonomy(
            logger,
            args.contig_fasta,
            args.cannot_name,
            args.GTDB_reference,
            num_cpu,
            binned_length,
            must_link_threshold,
            output)
        logger.info('Training model and clustering.')
        training(logger, [args.contig_fasta], args.bams,
                 args.num_process, [data_path], [data_split_path],
                 [os.path.join(output, 'cannot', 'cannot.txt')],
                 args.batchsize, args.epoches, output, device, args.ratio, args.min_len,  mode='single')

    binning(logger, args.bams, args.num_process, data_path,
            args.max_edges, args.max_node, args.minfasta_kb * 1000,
            binned_length, contig_length_dict, contig_dict,recluster,
            os.path.join(output, 'model.h5') if environment is None else None,
            random_seed, output,  device, environment if environment is not None else None)


def multi_easy_binning(args, logger, recluster,
                       random_seed, output, device):
    """
    contain `predict_taxonomy`, `generate_data_multi`, `train`, `bin` in one command for multi-sample binning
    """
    logger.info('Multi-sample binning.')
    logger.info('Generate training data.')

    num_cpu = multiprocessing.cpu_count() if args.num_process == 0 else args.num_process

    sample_list = generate_data_multi(
        logger,
        args.contig_fasta,
        args.bams,
        args.num_process,
        args.separator,
        args.ratio,
        args.min_len,
        args.ml_threshold,
        output, )

    for sample in sample_list:
        logger.info(
            'Running mmseqs and generate cannot-link file of {}.'.format(sample))
        sample_fasta = os.path.join(
            output, 'samples', '{}.fasta'.format(sample))
        sample_data = os.path.join(output, 'samples', sample, 'data.csv')
        sample_data_split = os.path.join(
            output, 'samples', sample, 'data_split.csv')

        binned_short, must_link_threshold , contig_length_dict, contig_dict = process_fasta(sample_fasta, args.ratio)

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
            args.GTDB_reference,
            num_cpu,
            binned_length,
            must_link_threshold,
            os.path.join(output, 'samples', sample))

        sample_cannot = os.path.join(
            output, 'samples', sample, 'cannot/{}.txt'.format(sample))
        logger.info('Training model and clustering for {}.'.format(sample))
        training(logger, [sample_fasta], args.bams, args.num_process,
                 [sample_data], [sample_data_split], [sample_cannot],
                 args.batchsize, args.epoches, os.path.join(output, 'samples', sample),
                 device, args.ratio, args.min_len, mode='single')

        binning(logger, args.bams, args.num_process, sample_data,
                args.max_edges, args.max_node, args.minfasta_kb * 1000,
                binned_length, contig_length_dict, contig_dict, recluster,
                os.path.join(output, 'samples', sample, 'model.h5'), random_seed ,
                os.path.join(output, 'samples', sample),  device, None)

    os.makedirs(os.path.join(output, 'bins'), exist_ok=True)
    for sample in sample_list:
        bin_dir_name = 'output_recluster_bins' if recluster else 'output_bins'
        for bf in os.listdir(os.path.join(
            output, 'samples', sample, bin_dir_name)):
            original_path = os.path.join(
                output, 'samples', sample, bin_dir_name, bf)
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

    validate_args(args)

    if args.cmd != 'download_GTDB':
        import torch
        out = args.output
        os.makedirs(out, exist_ok=True)

        device = torch.device(
            "cuda" if torch.cuda.is_available() else "cpu")

        if args.cmd != 'train':
            if os.path.splitext(args.contig_fasta)[1] == '.gz':
                contig_name = unzip_fasta('gz', args.contig_fasta)
                args.contig_fasta = contig_name
            elif os.path.splitext(args.contig_fasta)[1] == '.bz2':
                contig_name = unzip_fasta('bz2', args.contig_fasta)
                args.contig_fasta = contig_name
        else:
            contig_fastas = []
            for contig in args.contig_fasta:
                if os.path.splitext(contig)[1] == '.gz':
                    contig_name = unzip_fasta('gz', args.contig_fasta)
                    contig_fastas.append(contig_name)
                elif os.path.splitext(contig)[1] == '.bz2':
                    contig_name = unzip_fasta('bz2', args.contig_fasta)
                    contig_fastas.append(contig_name)
                else:
                    contig_fastas.append(contig)
            args.contig_fasta = contig_fastas

    if args.cmd in ['predict_taxonomy', 'generate_data_single', 'bin','single_easy_bin']:
        binned_short, must_link_threshold, contig_length_dict, contig_dict = process_fasta(args.contig_fasta, args.ratio)
        if args.min_len is None:
            binned_length = 1000 if binned_short else 2500
        else:
            binned_length = args.min_len

    if args.cmd in ['predict_taxonomy', 'generate_data_single', 'generate_data_multi', 'single_easy_bin', 'multi_easy_bin']:
        if args.ml_threshold is not None:
            must_link_threshold = args.ml_threshold

    if args.cmd == 'download_GTDB':
        download_GTDB(logger,args.GTDB_reference)

    if args.cmd == 'predict_taxonomy':
        predict_taxonomy(
            logger,
            args.contig_fasta,
            args.cannot_name,
            args.GTDB_reference,
            args.num_process,
            binned_length,
            must_link_threshold,
            out)

    if args.cmd == 'generate_data_single':
        generate_data_single(
            logger,
            args.contig_fasta,
            args.bams,
            binned_length,
            must_link_threshold,
            args.num_process,
            out)

    if args.cmd == 'generate_data_multi':
        generate_data_multi(
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
        training(logger, args.contig_fasta, args.bams, args.num_process,
                 args.data, args.data_split, args.cannot_link,
                 args.batchsize, args.epoches, out, device, args.ratio, args.min_len, args.mode)


    if args.cmd == 'bin':
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        binning(logger,args.bams, args.num_process, args.data, args.max_edges,
                args.max_node, args.minfasta_kb * 1000, binned_length,contig_length_dict,
                contig_dict, args.recluster, args.model_path, args.random_seed,out, device, args.environment)


    if args.cmd == 'single_easy_bin':
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        single_easy_binning(
            args,
            logger,
            binned_length,
            must_link_threshold,
            contig_length_dict,
            contig_dict,
            args.recluster,
            args.random_seed,
            out,
            device, args.environment)

    if args.cmd == 'multi_easy_bin':
        if args.random_seed is not None:
            set_random_seed(args.random_seed)
        multi_easy_binning(
            args,
            logger,
            args.recluster,
            args.random_seed,
            out,
            device)


if __name__ == '__main__':
    main()
