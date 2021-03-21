import argparse
import sys
import logging
import os
from .utils import validate_args, get_threshold, generate_cannot_link
import subprocess
import gzip
import bz2
from Bio import SeqIO
import multiprocessing
from .generate_coverage import calculate_coverage
from atomicwrites import atomic_write
from .generate_kmer import generate_kmer_features_from_fasta
import pandas as pd
from Bio.SeqRecord import SeqRecord
from .semi_supervised_model import train
import torch
from .cluster import cluster
import shutil

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural network for metagenomic binning')

    subparsers = parser.add_subparsers(title='S3N2Bin subcommands',
                                       dest='cmd',
                                       metavar='')

    single_easy_bin = subparsers.add_parser('single_easy_bin',
                                            help='Get the bins with single or co-assembly binning using one line command.')

    multi_easy_bin = subparsers.add_parser('multi_easy_bin',
                                            help='Get the bins with multi-samples binning using one line command.')

    predict_taxonomy = subparsers.add_parser('predict_taxonomy',
                                             help='Run the contig annotation using mmseqs '
                                                  'with GTDB reference genome and generate'
                                                  'cannot-link file used in the semi-supervsied deep learning model training.'
                                                  '(Will download the GTDB database if not input path of GTDB )')

    generate_data_single = subparsers.add_parser('generate_data_single',
                                                 help='Generate training data(data.csv,data_split.csv) '
                                                      'for single and co-assembly binning '
                                                      'for the semi-supervised deep learning model training.')

    generate_data_multi = subparsers.add_parser('generate_data_multi', help='Generate training data(data.csv,data_split.csv) '
                                          'for multi-samples binning '
                                          'for the semi-supervised deep learning model training.')

    binning = subparsers.add_parser('bin',
                                    help='Training the model and clustering contigs to bins.')

    binning.add_argument('--data',
                         required=True,
                         help='Path to the input data.csv file.',
                         dest='data',
                         default=None,
                         )
    binning.add_argument('--data-split',
                         required=True,
                         help='Path to the input data_split.csv file.',
                         dest='data_split',
                         default=None,
                         )
    binning.add_argument('-c', '--cannot-link',
                         required=True,
                         nargs='*',
                         help='Path to the input cannot link file generated from other additional biological information,'
                         'one row for each cannot link constraint.'
                         'The file format: contig_1,contig_2.',
                         dest='cannot_link',
                         default=None,
                         metavar='')

    for p in [single_easy_bin, multi_easy_bin, predict_taxonomy, generate_data_single, generate_data_multi,binning]:
        p.add_argument('-i', '--input-fasta',
                                required=True,
                                help='Path to the input fasta file.',
                                dest='contig_fasta',
                                default=None,)

    for p in [single_easy_bin, multi_easy_bin, predict_taxonomy, generate_data_single, generate_data_multi, binning]:
        p.add_argument('-o', '--output',
                            required=True,
                            help='Output directory (will be created if non-existent)',
                            dest='output',
                            default=None,
                            )

    for p in [single_easy_bin, multi_easy_bin, generate_data_single, generate_data_multi, binning]:
        p.add_argument('-b', '--input-bam',
                            required=True,
                            nargs='*',
                            help='Path to the input BAM file. '
                                 'If using multiple sample binning, you can input multiple files.',
                            dest='bams',
                            default=None,
                            )

        p.add_argument('-p', '--processes', '-t', '--threads',
                                     required=False,
                                     type=int,
                                     help='Number of CPUs used (pass the value 0 to use all CPUs)',
                                     dest='num_process',
                                     default=0,
                                     metavar=''
                                     )

    for p in [single_easy_bin, multi_easy_bin, predict_taxonomy]:
        p.add_argument('-r', '--reference-db',
                            required=False,
                            help='GTDB reference file. (Default: $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB).'
                            'If not set --reference-db and can not find GTDB in $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB, '
                            'we will download GTDB to the default path.',
                            dest='GTDB_reference',
                            metavar='',
                            default=None)
    for p in [single_easy_bin, predict_taxonomy]:
        p.add_argument('--cannot-name',
                            required=False,
                            help='Name for the cannot-link file.',
                            dest='cannot_name',
                            default='cannot',
                            metavar=''
                            )

    for p in [single_easy_bin, multi_easy_bin, binning]:
        p.add_argument('--minfasta-kbs',
                            required=False,
                            type=int,
                            help='minimum bin size in Kbps (Default: 200).',
                            dest='minfasta_kb',
                            default=200,
                            metavar='')

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


    if not args:
        parser.print_help(sys.stderr)
        sys.exit()
    return parser.parse_args(args)


def generate_cov(bam_file, bam_index, out, threshold,
                 is_combined, contig_threshold, logger):
    logger.info('Processing {}'.format(bam_file))
    bam_name = os.path.split(bam_file)[-1] + '_{}'.format(bam_index)
    bam_depth = os.path.join(out, '{}_depth.txt'.format(bam_name))
    os.system(
        'bedtools genomecov -bga -ibam {0} > {1}'.format(bam_file, bam_depth))

    if is_combined:
        contig_cov, must_link_contig_cov = calculate_coverage(bam_depth, threshold, is_combined=is_combined,
                                                              contig_threshold=contig_threshold)
        contig_cov = contig_cov.apply(lambda x: x + 1e-5)
        must_link_contig_cov = must_link_contig_cov.apply(lambda x: x + 1e-5)
        contig_cov = contig_cov / 100
        must_link_contig_cov = must_link_contig_cov / 100
        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)

        with atomic_write(os.path.join(out, '{}_data_split_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            must_link_contig_cov.to_csv(ofile)
    else:
        contig_cov = calculate_coverage(
            bam_depth,
            threshold,
            is_combined=is_combined,
            contig_threshold=contig_threshold)
        contig_cov = contig_cov.apply(lambda x: x + 1e-5)
        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)
    return (bam_file, logger)


def generate_cov_multiple(bam_file, bam_index, out, threshold,
                          is_combined, sep, binned_threshold_dict, logger):
    logger.info('Processing {}'.format(bam_file))
    bam_name = os.path.split(bam_file)[-1] + '_{}'.format(bam_index)
    bam_depth = os.path.join(out, '{}_depth.txt'.format(bam_name))
    os.system(
        'bedtools genomecov -bga -ibam {0} > {1}'.format(bam_file, bam_depth))
    if is_combined:
        contig_cov, must_link_contig_cov = calculate_coverage(bam_depth, threshold, is_combined=is_combined,
                                                              sep=sep, binned_thre_dict=binned_threshold_dict)
        contig_cov = contig_cov.apply(lambda x: x + 1e-5)
        must_link_contig_cov = must_link_contig_cov.apply(lambda x: x + 1e-5)
        contig_cov = contig_cov / 100
        must_link_contig_cov = must_link_contig_cov / 100
        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)

        with atomic_write(os.path.join(out, '{}_data_split_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            must_link_contig_cov.to_csv(ofile)
    else:
        contig_cov = calculate_coverage(bam_depth, threshold, is_combined=is_combined, sep=sep,
                                        binned_thre_dict=binned_threshold_dict)
        contig_cov = contig_cov.apply(lambda x: x + 1e-5)
        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)
    return (bam_file, logger)


def _checkback(msg):
    msg[1].info('Processed:{}'.format(msg[0]))


def predict_taxonomy(contig_fasta, GTDB_reference,
                     cannot_name, logger,
                     output, handle,
                     binned_short, must_link_threshold):
    """
    Predict taxonomy using mmseqs and generate cannot-link file
    :param handle: handle to read fasta file
    :param binned_short: binning contigs > 1000bp or 2500 bp
    :param must_link_threshold: threshold for generating must-link pair
    """
    GTDB_default = os.path.join(
        os.environ['HOME'],
        '.cache',
        'S3N2Bin',
        'mmseqs2-GTDB',
        'GTDB')
    GTDB_path = GTDB_reference
    if GTDB_reference is None:
        if not os.path.exists(GTDB_default):
            logger.info('Downloading GTDB.')
            GTDB_dir = os.path.split(GTDB_default)[0]
            os.makedirs(GTDB_dir, exist_ok=True)
            subprocess.check_call(
                ['mmseqs',
                 'databases',
                 'GTDB',
                 GTDB_default,
                 '{}/tmp'.format(GTDB_dir),
                 ],
                stdout=None,
                stderr=subprocess.DEVNULL,
            )
        GTDB_path = GTDB_default
    subprocess.check_call(
        ['mmseqs',
         'createdb',
         contig_fasta,
         '{}/contig_DB'.format(output)],
        stdout=None,
        stderr=subprocess.DEVNULL,
    )
    os.makedirs(os.path.join(output, 'mmseqs_annotation'), exist_ok=True)
    subprocess.run(
        ['mmseqs',
         'taxonomy',
         '{}/contig_DB'.format(output),
         GTDB_path,
         os.path.join(output, 'mmseqs_annotation/mmseqs_annotation'),
         os.path.join(output, 'mmseqs_tmp'),
         '--tax-lineage',
         str(1),
         ],
        check=True,
        stdout=None,
        stderr=subprocess.DEVNULL,
    )
    subprocess.check_call(
        ['mmseqs',
         'createtsv',
         '{}/contig_DB'.format(output),
         os.path.join(output, 'mmseqs_annotation/mmseqs_annotation'),
         os.path.join(output, 'mmseqs_annotation/taxonomyResult.tsv')
         ],
        stdout=None,
        stderr=subprocess.DEVNULL,
    )

    namelist = []
    num_must_link = 0
    binned_threshold = 1000 if binned_short else 2500
    for seq_record in SeqIO.parse(handle, "fasta"):
        if len(seq_record) > binned_threshold:
            namelist.append(seq_record.id)
        if len(seq_record) >= must_link_threshold:
            num_must_link += 1
    os.makedirs(os.path.join(output, 'cannot'), exist_ok=True)
    generate_cannot_link(
        os.path.join(output, 'mmseqs_annotation/taxonomyResult.tsv'),
        namelist, num_must_link,
        os.path.join(output, 'cannot'), cannot_name)


def generate_data_single(bams, num_process, logger,
                         output, handle, binned_short,
                         must_link_threshold):
    """
    Generate data.csv and data_split.csv for training and clustering.
    data.csv has the features(kmer and abundance) for original contigs.
    data_split.csv has the features(kmer and abundace) for contigs that are breaked up as must-link pair.
    """
    n_sample = len(bams)
    is_combined = n_sample >= 5
    bam_list = bams
    if num_process != 0:
        pool = multiprocessing.Pool(num_process)
    else:
        pool = multiprocessing.Pool()

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
                1000 if binned_short else 2500,
                logger,
            ),
            callback=_checkback)
    pool.close()
    pool.join()

    logger.info('Start generating kmer features from fasta file.')
    kmer_whole = generate_kmer_features_from_fasta(
        handle, 1000 if binned_short else 2500, 4)
    kmer_split = generate_kmer_features_from_fasta(
        handle, 1000, 4, split=True, threshold=must_link_threshold)

    data = kmer_whole
    data_split = kmer_split
    data.index = data.index.astype(str)

    for bam_index, bam_file in enumerate(bam_list):
        cov = pd.read_csv(os.path.join(output, '{}_data_cov.csv'.format(
            os.path.split(bam_file)[-1] + '_{}'.format(bam_index))), index_col=0)
        cov.index = cov.index.astype(str)
        data = pd.merge(data, cov, how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)
        if is_combined:
            cov_split = pd.read_csv(os.path.join(output, '{}_data_split_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index))), index_col=0)

            data_split = pd.merge(data_split, cov_split, how='inner', on=None,
                                  left_index=True, right_index=True, sort=False, copy=True)
    if not is_combined:
        data_split = kmer_split

    with atomic_write(os.path.join(output, 'data.csv'), overwrite=True) as ofile:
        data.to_csv(ofile)

    with atomic_write(os.path.join(output, 'data_split.csv'), overwrite=True) as ofile:
        data_split.to_csv(ofile)


def generate_data_multi(bams, num_process,separator,
                        logger, output, handle):
    n_sample = len(bams)
    is_combined = n_sample >= 5
    bam_list = bams
    if num_process != 0:
        pool = multiprocessing.Pool(num_process)
    else:
        pool = multiprocessing.Pool()

    # Gererate contig file for every sample
    from collections import defaultdict
    sample_list = list()
    contig_sample_list = []
    contig_length_list = []
    flag_name = None

    os.makedirs(os.path.join(output, 'samples'), exist_ok=True)

    for seq_record in SeqIO.parse(handle, "fasta"):
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

    must_link_threshold = get_threshold(contig_length_list)
    logger.info('Calculating coverage for every sample.')

    binning_threshold = {1000: [], 2500: []}
    for sample in sample_list:
        whole_contig_bp = 0
        contig_bp_2500 = 0
        for seq_record in SeqIO.parse(os.path.join(
                output, 'samples/{}.fasta'.format(sample)), "fasta"):
            if len(seq_record) >= 1000 and len(seq_record) <= 2500:
                contig_bp_2500 += len(seq_record)
            whole_contig_bp += len(seq_record)
        binned_short = contig_bp_2500 / whole_contig_bp < 0.05
        binning_threshold[1000].append(sample) if binned_short \
            else binning_threshold[2500].append(sample)

    for bam_index in range(n_sample):
        pool.apply_async(generate_cov_multiple,
                         args=(
                             bam_list[bam_index],
                             bam_index,
                             os.path.join(output, 'samples'),
                             must_link_threshold,
                             is_combined,
                             separator,
                             binning_threshold,
                             logger
                         ),
                         callback=_checkback)
    pool.close()
    pool.join()
    # Generate cov features for every sample
    data_cov = pd.read_csv(os.path.join(output, 'samples', '{}_data_cov.csv'.format(
        os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)
    if is_combined:
        data_split_cov = pd.read_csv(os.path.join(
            output, 'samples', '{}_data_split_cov.csv'.format(
                os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)
    data_cov.index = data_cov.index.astype(str)
    for bam_index, bam_file in enumerate(bam_list):
        if bam_index == 0:
            continue
        cov = pd.read_csv(
            os.path.join(output, 'samples', '{}_data_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index))),
            index_col=0)
        cov.index = cov.index.astype(str)
        data_cov = pd.merge(data_cov, cov, how='inner', on=None,
                            left_index=True, right_index=True, sort=False, copy=True)

        if is_combined:
            cov_split = pd.read_csv(os.path.join(output, 'samples', '{}_data_split_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index))), index_col=0)

            data_split_cov = pd.merge(data_split_cov, cov_split, how='inner', on=None,
                                      left_index=True, right_index=True, sort=False, copy=True)

    data_cov = data_cov.reset_index()
    columns_list = list(data_cov.columns)
    columns_list[0] = 'contig_name'
    data_cov.columns = columns_list

    if is_combined:
        data_split_cov = data_split_cov.reset_index()
        columns_list = list(data_split_cov.columns)
        columns_list[0] = 'contig_name'
        data_split_cov.columns = columns_list

    for sample in sample_list:
        output_path = os.path.join(output, 'samples', sample)
        os.makedirs(output_path, exist_ok=True)
        part_data = data_cov[data_cov['contig_name'].str.contains(
            '{}'.format(sample + separator))]
        part_data = part_data.set_index('contig_name')
        part_data.index.name = None
        index_list = part_data.index.tolist()
        index_list = [temp.split(separator)[1] for temp in index_list]
        part_data.index = index_list
        part_data.to_csv(os.path.join(output_path, 'data_cov.csv'))
        if is_combined:
            part_data = data_split_cov[data_split_cov['contig_name'].str.contains(
                '{}'.format(sample + separator))]
            part_data = part_data.set_index('contig_name')
            part_data.index.name = None
            index_list = part_data.index.tolist()
            index_list = [temp.split(separator)[1]
                          for temp in index_list]
            part_data.index = index_list
            part_data.to_csv(os.path.join(
                output_path, 'data_split_cov.csv'))

        sample_contig_fasta = os.path.join(
            output, 'samples/{}.fasta'.format(sample))
        binned_short = True if sample in binning_threshold[1000] else False
        kmer_whole = generate_kmer_features_from_fasta(
            sample_contig_fasta, 1000 if binned_short else 2500, 4)
        kmer_split = generate_kmer_features_from_fasta(
            sample_contig_fasta, 1000, 4, split=True, threshold=must_link_threshold)


        sample_cov = pd.read_csv(os.path.join(
            output_path, 'data_cov.csv'),
            index_col=0)
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


def binning(contig_fasta, bams, num_process, data,
            data_split, cannot_link, batchsize, epoches,
            max_edges, max_node, minfasta, logger, output,
            binned_short, device, contig_length_dict, contig_dict):
    """
    Training and clustering the contigs to get the final bins.
    """
    logger.info('Start binning.')
    n_sample = len(bams)
    is_combined = n_sample >= 5
    num_cpu = multiprocessing.cpu_count() if num_process == 0 else num_process
    data = pd.read_csv(data, index_col=0)
    data.index = data.index.astype(str)
    data_split = pd.read_csv(data_split, index_col=0)

    model = train(
        output,
        contig_fasta,
        binned_short,
        logger,
        data,
        data_split,
        cannot_link,
        is_combined,
        batchsize,
        epoches,
        device,
        num_cpu)
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
        binned_short,
        num_cpu,
        minfasta)


def single_easy_binning(args, logger, output, handle, binned_short,
                        must_link_threshold, device, contig_length_dict, contig_dict):
    logger.info('Running mmseqs and generate cannot-link file.')
    predict_taxonomy(
        args.contig_fasta,
        args.GTDB_reference,
        args.cannot_name,
        logger,
        output, handle, binned_short, must_link_threshold)
    logger.info('Generate training data.')
    generate_data_single(
        args.bams,
        args.num_process,
        logger,
        output, handle, binned_short, must_link_threshold)
    logger.info('Training model and clustering.')
    data_path = os.path.join(output,'data.csv')
    data_split_path = os.path.join(output,'data_split.csv')
    binning(args.contig_fasta, args.bams, args.num_process, data_path,
            data_split_path, os.path.join(
                output, 'cannot', 'cannot.txt'), args.batchsize, args.epoches,
            args.max_edges, args.max_node, args.minfasta_kb * 1000, logger, output, binned_short, device, contig_length_dict, contig_dict)


def multi_easy_binning(args, logger, output, handle, device):
    logger.info('Multi-samples binning.')
    logger.info('Generate training data.')
    sample_list = generate_data_multi(
        args.bams,
        args.num_process,
        args.separator,
        logger,
        output, handle)

    for sample in sample_list:
        logger.info(
            'Running mmseqs and generate cannot-link file of {}.'.format(sample))
        sample_fasta = os.path.join(
            output, 'samples', '{}.fasta'.format(sample))
        sample_data = os.path.join(output, 'samples', sample, 'data.csv')
        sample_data_split = os.path.join(
            output, 'samples', sample, 'data_split.csv')

        whole_contig_bp = 0
        contig_bp_2500 = 0
        contig_length_list = []
        contig_length_dict = {}
        contig_dict = {}
        handle = sample_fasta
        for seq_record in SeqIO.parse(handle, "fasta"):
            if len(seq_record) >= 1000 and len(seq_record) <= 2500:
                contig_bp_2500 += len(seq_record)
            contig_length_list.append(len(seq_record))
            whole_contig_bp += len(seq_record)
            contig_length_dict[str(seq_record.id).strip(
                '')] = len((seq_record.seq))
            contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

        binned_short = contig_bp_2500 / whole_contig_bp < 0.05
        must_link_threshold = get_threshold(contig_length_list)

        predict_taxonomy(
            sample_fasta,
            args.GTDB_reference,
            sample,
            logger,
            os.path.join(output, 'samples', sample), handle, binned_short, must_link_threshold,)
        sample_cannot = os.path.join(
            output, 'samples', sample, 'cannot/{}.txt'.format(sample))
        logger.info('Training model and clustering for {}.'.format(sample))
        binning(sample_fasta, args.bams, args.num_process, sample_data,
                sample_data_split, sample_cannot, args.batchsize, args.epoches,
                args.max_edges, args.max_node, args.minfasta_kb * 1000, logger, os.path.join(output, 'samples', sample), binned_short, device, contig_length_dict, contig_dict)

    os.makedirs(os.path.join(output, 'bins'), exist_ok=True)
    for sample in sample_list:
        bin_file = os.listdir(os.path.join(
            output, 'samples', sample, 'output_recluster_bins'))
        for bin in bin_file:
            original_path = os.path.join(
                output, 'samples', sample, 'output_recluster_bins', bin)
            new_file = '{0}_{1}'.format(sample, bin)
            new_path = os.path.join(output, 'bins', new_file)
            shutil.copyfile(original_path, new_path)


def main():
    args = sys.argv[1:]
    args = parse_args(args)

    logger = logging.getLogger('S3N2Bin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    validate_args(args)

    out = args.output
    os.makedirs(out, exist_ok=True)

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu")

    if os.path.splitext(args.contig_fasta)[1] == '.gz':
        handle = gzip.open(args.contig_fasta, "rt")
    elif os.path.splitext(args.contig_fasta)[1] == '.bz2':
        handle = bz2.open(args.contig_fasta, "rt")
    else:
        handle = args.contig_fasta

    if args.cmd in ['predict_taxonomy', 'generate_data_single', 'bin','single_easy_bin']:
        whole_contig_bp = 0
        contig_bp_2500 = 0
        contig_length_list = []
        contig_length_dict = {}
        contig_dict = {}

        for seq_record in SeqIO.parse(handle, "fasta"):
            if len(seq_record) >= 1000 and len(seq_record) <= 2500:
                contig_bp_2500 += len(seq_record)
            contig_length_list.append(len(seq_record))
            whole_contig_bp += len(seq_record)
            contig_length_dict[str(seq_record.id).strip(
                '')] = len((seq_record.seq))
            contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

        binned_short = contig_bp_2500 / whole_contig_bp < 0.05
        must_link_threshold = get_threshold(contig_length_list)

    if args.cmd == 'predict_taxonomy':
        predict_taxonomy(
            args.contig_fasta,
            args.GTDB_reference,
            args.cannot_name,
            logger,
            out,
            handle, binned_short, must_link_threshold)

    if args.cmd == 'generate_data_single':
        generate_data_single(
            args.bams,
            args.num_process,
            logger,
            out,
            handle, binned_short, must_link_threshold)

    if args.cmd == 'generate_data_multi':
        generate_data_multi(
            args.bams,
            args.num_process,
            args.separator,
            logger,
            out,
            handle)

    if args.cmd == 'bin':
        binning(args.contig_fasta, args.bams,
                args.num_process, args.data,
                args.data_split, args.cannot_link[0],
                args.batchsize, args.epoches,
                args.max_edges, args.max_node,
                args.minfasta_kb * 1000, logger,
                out, binned_short, device,
                contig_length_dict, contig_dict)

    if args.cmd == 'single_easy_bin':
        single_easy_binning(
            args,
            logger,
            out,
            handle,
            binned_short,
            must_link_threshold,
            device,
            contig_length_dict,
            contig_dict)

    if args.cmd == 'multi_easy_bin':
        multi_easy_binning(
            args,
            logger,
            out,
            handle,
            device)


if __name__ == '__main__':
    main()
