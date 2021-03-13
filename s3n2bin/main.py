import argparse
import os
import sys
import logging
from .generate_kmer import generate_kmer_features_from_fasta
from .generate_coverage import calculate_coverage
from Bio import SeqIO
from .utils import validate_args, get_threshold, generate_cannot_link
import pandas as pd
from .semi_supervised_model import train
import torch
from atomicwrites import atomic_write
from .cluster import cluster
import warnings
import multiprocessing
from Bio.SeqRecord import SeqRecord
import gzip
import bz2
import shutil
import subprocess


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural network for metagenomic binning')

    basic = parser.add_argument_group(title='Basic commands', description=None)
    basic.add_argument('command',
                       nargs=1,
                       help='You can choose easy-bin mode to get results with one line code (for single and co-assembly binning.)')

    basic.add_argument('-i', '--input-fasta',
                       required=True,
                       help='Path to the input contig fasta file.',
                       dest='contig_fasta',
                       default=None)
    basic.add_argument('-b', '--input-bam',
                       required=True,
                       nargs='*',
                       help='Path to the input bam file. '
                             'If using multiple sample binning, you can input multiple files.',
                       dest='bams',
                       default=None)
    basic.add_argument('-c', '--cannot-link',
                       required=False,
                       nargs='*',
                       help='Path to the input cannot link file generated from other additional biological information,'
                             'one row for each cannot link constraint.'
                             'The file format: contig_1,contig_2.',
                       dest='cannot_link',
                       default=None)
    basic.add_argument('-o', '--output',
                       required=True,
                       help='Output directory (will be created if non-existent)',
                       dest='output',
                       default=None)
    basic.add_argument('-s', '--separator',
                       required=False,
                       type=str,
                       help='Used when multiple samples binning to separate sample name and contig name.'
                            '(None means single sample and coassemble binning)',
                       dest='separator',
                       default=None,
                       )
    basic.add_argument('--GTDB-path',
                       required=False,
                       help='Path to the GTDB database used in the mmseqs annotation in easy-bin mode.'
                            '(If not set, we will download GTDB dataset to the output folder)',
                       dest='gtdb_path',
                       default=None,
                       )

    optional = parser.add_argument_group(
        title='Optional commands', description=None)
    optional.add_argument('-p', '--processes', '-t', '--threads',
                          required=False,
                          type=int,
                          help='Number of CPUs used',
                          dest='num_process',
                          default=0)

    optional.add_argument('--epoches',
                          required=False,
                          type=int,
                          help='Number of epoches used in the training process.',
                          dest='epoches',
                          default=20)

    optional.add_argument('--batch-size',
                          required=False,
                          type=int,
                          help='Batch size used in the training process.',
                          dest='batchsize',
                          default=2048,)

    optional.add_argument('--max-edges',
                          required=False,
                          type=int,
                          help='The maximum number of edges that can be connected to one contig.',
                          dest='max_edges',
                          default=200)

    optional.add_argument('--max-node',
                          required=False,
                          type=float,
                          dest='max_node',
                          default=1,
                          help='Percentage of contigs that considered to be binned.')

    optional.add_argument('--generate-data',
                          help='Used when multi-samples binning.'
                               'S3N2Bin will only generate datas (data.csv,data_split.csv) for training and clustering.',
                          required=False,
                          action='store_true',
                          dest='generate_data'
                          )
    optional.add_argument('--split-run',
                          help='Used when multi-samples binning after generating datas for training and clustering.'
                               'With this command you can run S3N2Bin parallelly on samples manually.',
                          required=False,
                          action='store_true',
                          dest='split_running'
                          )

    return parser.parse_args()


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


def main(args=None):
    if args is None:
        args = sys.argv

    args = parse_args(args)
    [args.command] = args.command
    validate_args(args)
    logger = logging.getLogger('S3N2Bin')
    logger.setLevel(logging.INFO)

    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    out = args.output
    os.makedirs(out, exist_ok=True)

    contig_fasta = args.contig_fasta

    multi_binning = False
    if args.separator is not None:
        multi_binning = True

    if args.split_running:
        multi_binning = False

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu")

    if not multi_binning:
        whole_contig_bp = 0
        contig_bp_2500 = 0
        contig_length_list = []
        contig_length_dict = {}
        contig_dict = {}
        if os.path.splitext(contig_fasta)[1] == '.gz':
            handle = gzip.open(contig_fasta, "rt")
        elif os.path.splitext(contig_fasta)[1] == '.bz2':
            handle = bz2.open(contig_fasta, "rt")
        else:
            handle = contig_fasta
        for seq_record in SeqIO.parse(handle, "fasta"):
            if len(seq_record) >= 1000 and len(seq_record) <= 2500:
                contig_bp_2500 += len(seq_record)
            contig_length_list.append(len(seq_record))
            whole_contig_bp += len(seq_record)
            contig_length_dict[str(seq_record.id).strip(
                '')] = len((seq_record.seq))
            contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

        binned_short = contig_bp_2500 / whole_contig_bp < 0.05
        n_sample = len(args.bams)
        is_combined = n_sample >= 5

        # threshold for generating must link pairs
        threshold = get_threshold(contig_length_list)

        if args.command == 'easy-bin':
            logger.info('Running mmseqs annotation')
            if args.gtdb_path is None:
                logger.info('Downloading GTDB.')
                subprocess.check_call(
                    ['mmseqs',
                     'databases',
                     'GTDB',
                     '{}/GTDB'.format(out),
                     '{}/tmp'.format(out),
                     ],
                    stdout=None,
                    stderr=subprocess.DEVNULL,
                )
            gtdb_path = os.path.join(
                out, 'GTDB') if args.gtdb_path is None else args.gtdb_path
            subprocess.check_call(
                ['mmseqs',
                 'createdb',
                 contig_fasta,
                 '{}/contig_DB'.format(out)],
                stdout=None,
                stderr=subprocess.DEVNULL,
            )
            os.makedirs(os.path.join(out, 'mmseqs_annotation'), exist_ok=True)
            subprocess.run(
                ['mmseqs',
                 'taxonomy',
                 '{}/contig_DB'.format(out),
                 gtdb_path,
                 os.path.join(out, 'mmseqs_annotation/mmseqs_annotation'),
                 os.path.join(out, 'mmseqs_tmp'),
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
                 '{}/contig_DB'.format(out),
                 os.path.join(out, 'mmseqs_annotation/mmseqs_annotation'),
                 os.path.join(out, 'mmseqs_annotation/taxonomyResult.tsv')
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
                if len(seq_record) >= threshold:
                    num_must_link += 1
            os.makedirs(os.path.join(out, 'cannot'), exist_ok=True)
            generate_cannot_link(
                os.path.join(
                    out,
                    'mmseqs_annotation/taxonomyResult.tsv'),
                namelist,
                num_must_link,
                os.path.join(
                    out,
                    'cannot'),
                'cannot')

        if not args.split_running:
            logger.info('Calculating coverage for every sample.')
            # generating coverage for every contig and for must link pair

            if args.num_process != 0:
                pool = multiprocessing.Pool(args.num_process)
            else:
                pool = multiprocessing.Pool()

            bam_list = args.bams
            for bam_index in range(n_sample):
                pool.apply_async(
                    generate_cov,
                    args=(
                        bam_list[bam_index],
                        bam_index,
                        out,
                        threshold,
                        is_combined,
                        1000 if binned_short else 2500,
                        logger,
                    ),
                    callback=_checkback)
            pool.close()
            pool.join()

            # Processing input contig fasta file
            logger.info('Start generating kmer features from fasta file.')
            kmer_whole = generate_kmer_features_from_fasta(
                args.contig_fasta, 1000 if binned_short else 2500, 4)
            kmer_split = generate_kmer_features_from_fasta(
                args.contig_fasta, 1000, 4, split=True, threshold=threshold)

            data = kmer_whole
            data_split = kmer_split

            for bam_index, bam_file in enumerate(bam_list):
                cov = pd.read_csv(os.path.join(out, '{}_data_cov.csv'.format(
                    os.path.split(bam_file)[-1] + '_{}'.format(bam_index))), index_col=0)
                data = pd.merge(data, cov, how='inner', on=None,
                                left_index=True, right_index=True, sort=False, copy=True)
                if is_combined:
                    cov_split = pd.read_csv(os.path.join(out, '{}_data_split_cov.csv'.format(
                        os.path.split(bam_file)[-1] + '_{}'.format(bam_index))), index_col=0)

                    data_split = pd.merge(data_split, cov_split, how='inner', on=None,
                                          left_index=True, right_index=True, sort=False, copy=True)
            if not is_combined:
                data_split = kmer_split

            with atomic_write(os.path.join(out, 'data.csv'), overwrite=True) as ofile:
                data.to_csv(ofile)

            with atomic_write(os.path.join(out, 'data_split.csv'), overwrite=True) as ofile:
                data_split.to_csv(ofile)
        else:
            data = pd.read_csv(os.path.join(out, 'data.csv'), index_col=0)
            data_split = pd.read_csv(
                os.path.join(out,'data_split.csv'),index_col=0)

        model = train(
            out,
            args.contig_fasta,
            binned_short,
            logger,
            data,
            data_split,
            args.cannot_link[0] if args.command != 'easy-bin' else os.path.join(
                out, 'cannot/cannot.txt'),
            is_combined,
            args.batchsize,
            args.epoches,
            device)
        cluster(
            model,
            data,
            device,
            args.max_edges,
            args.max_node,
            is_combined,
            logger,
            n_sample,
            contig_length_dict,
            out,
            contig_dict,
            binned_short)

    else:
        """
        Multi-samples binning
        """
        # Gererate contig file for every sample
        from collections import defaultdict
        sample_list = list()
        contig_sample_list = []
        contig_length_list = []
        flag_name = None

        os.makedirs(os.path.join(out, 'samples'), exist_ok=True)

        if os.path.splitext(contig_fasta)[1] == '.gz':
            handle = gzip.open(contig_fasta, "rt")
        elif os.path.splitext(contig_fasta)[1] == '.bz2':
            handle = bz2.open(contig_fasta, "rt")
        else:
            handle = contig_fasta

        for seq_record in SeqIO.parse(handle, "fasta"):
            sample_name, contig_name = seq_record.id.split(args.separator)
            if flag_name is None:
                flag_name = sample_name
            if sample_name == flag_name:
                rec = SeqRecord(seq_record.seq, id=contig_name, description='')
                contig_sample_list.append(rec)
            if sample_name != flag_name:
                SeqIO.write(contig_sample_list,
                            os.path.join(
                                out, 'samples', '{}.fasta'.format(flag_name)), 'fasta')
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
                            out, 'samples', '{}.fasta'.format(flag_name)), 'fasta')

        threshold = get_threshold(contig_length_list)
        logger.info('Calculating coverage for every sample.')

        binning_threshold = {1000: [], 2500: []}
        for sample in sample_list:
            whole_contig_bp = 0
            contig_bp_2500 = 0
            for seq_record in SeqIO.parse(os.path.join(
                    out, 'samples/{}.fasta'.format(sample)), "fasta"):
                if len(seq_record) >= 1000 and len(seq_record) <= 2500:
                    contig_bp_2500 += len(seq_record)
                whole_contig_bp += len(seq_record)
            binned_short = contig_bp_2500 / whole_contig_bp < 0.05
            binning_threshold[1000].append(sample) if binned_short \
                else binning_threshold[2500].append(sample)

        n_sample = len(args.bams)
        is_combined = n_sample >= 5
        bam_list = args.bams

        if args.num_process != 0:
            pool = multiprocessing.Pool(args.num_process)
        else:
            pool = multiprocessing.Pool()

        for bam_index in range(n_sample):
            pool.apply_async(generate_cov_multiple,
                             args=(
                                 bam_list[bam_index],
                                 bam_index,
                                 os.path.join(out, 'samples'),
                                 threshold,
                                 is_combined,
                                 args.separator,
                                 binning_threshold,
                                 logger
                             ),
                             callback=_checkback)
        pool.close()
        pool.join()

        # Generate cov features for every sample
        data_cov = pd.read_csv(os.path.join(out, 'samples', '{}_data_cov.csv'.format(
            os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)
        if is_combined:
            data_split_cov = pd.read_csv(os.path.join(
                out, 'samples', '{}_data_split_cov.csv'.format(
                    os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)
        for bam_index, bam_file in enumerate(bam_list):
            if bam_index == 0:
                continue
            cov = pd.read_csv(
                os.path.join(out, 'samples', '{}_data_cov.csv'.format(
                    os.path.split(bam_file)[-1] + '_{}'.format(bam_index))),
                index_col=0)
            data_cov = pd.merge(data_cov, cov, how='inner', on=None,
                                left_index=True, right_index=True, sort=False, copy=True)

            if is_combined:
                cov_split = pd.read_csv(os.path.join(out, 'samples', '{}_data_split_cov.csv'.format(
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
            output_path = os.path.join(out, 'samples', sample)
            os.makedirs(output_path, exist_ok=True)
            part_data = data_cov[data_cov['contig_name'].str.contains(
                '{}'.format(sample + args.separator))]
            part_data = part_data.set_index('contig_name')
            part_data.index.name = None
            index_list = part_data.index.tolist()
            index_list = [temp.split(args.separator)[1] for temp in index_list]
            part_data.index = index_list
            part_data.to_csv(os.path.join(output_path, 'data_cov.csv'))
            if is_combined:
                part_data = data_split_cov[data_split_cov['contig_name'].str.contains(
                    '{}'.format(sample + args.separator))]
                part_data = part_data.set_index('contig_name')
                part_data.index.name = None
                index_list = part_data.index.tolist()
                index_list = [temp.split(args.separator)[1]
                              for temp in index_list]
                part_data.index = index_list
                part_data.to_csv(os.path.join(
                    output_path, 'data_split_cov.csv'))

        cannot_link_dict = defaultdict(list)
        for cannot in args.cannot_link:
            cannot_file = os.path.split(cannot)[-1]
            cannot_link_dict[cannot_file.split('.')[0]] = cannot

        for sample in sample_list:
            logger.info('Clustering:{}'.format(sample))
            output_path = os.path.join(out, 'samples', sample)
            sample_contig_fasta = os.path.join(
                out, 'samples/{}.fasta'.format(sample))

            contig_length_dict = {}
            contig_dict = {}
            for seq_record in SeqIO.parse(sample_contig_fasta, "fasta"):
                contig_length_dict[str(seq_record.id).strip(
                    '')] = len((seq_record.seq))
                contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

            binned_short = True if sample in binning_threshold[1000] else False

            logger.info('Start generating kmer features from fasta file.')
            kmer_whole = generate_kmer_features_from_fasta(
                sample_contig_fasta, 1000 if binned_short else 2500, 4)
            kmer_split = generate_kmer_features_from_fasta(
                sample_contig_fasta, 1000, 4, split=True, threshold=threshold)

            sample_cov = pd.read_csv(os.path.join(
                output_path, 'data_cov.csv'),
                index_col=0)
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

            if not args.generate_data:
                model = train(
                    output_path,
                    sample_contig_fasta,
                    binned_short,
                    logger,
                    data,
                    data_split,
                    cannot_link_dict[sample],
                    is_combined,
                    args.batchsize,
                    args.epoches,
                    device)
                cluster(
                    model,
                    data,
                    device,
                    args.max_edges,
                    args.max_node,
                    is_combined,
                    logger,
                    n_sample,
                    contig_length_dict,
                    output_path,
                    contig_dict,
                    binned_short)

        if not args.generate_data:
            os.makedirs(os.path.join(out, 'bins'), exist_ok=True)
            for sample in sample_list:
                bin_file = os.listdir(os.path.join(
                    out, 'samples', sample, 'output_recluster_bins'))
                for bin in bin_file:
                    original_path = os.path.join(
                        out, 'samples', sample, 'output_recluster_bins', bin)
                    new_file = '{0}_{1}'.format(sample, bin)
                    new_path = os.path.join(out, 'bins', new_file)
                    shutil.copyfile(original_path, new_path)

    if __name__ == '__main__':
        warnings.filterwarnings('ignore')
        main(sys.argv)
