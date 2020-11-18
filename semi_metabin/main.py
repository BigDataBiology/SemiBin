import argparse
import os
import sys
import numpy as np
import logging
from .generate_kmer import generate_kmer_features_from_fasta
from .generate_coverage import calculate_coverage
from  Bio import SeqIO
import pandas as pd
import subprocess

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural network for metagenomic binning')
    parser.add_argument('-i','--input-fasta',
                        required=True,
                        help='Path to the input contig fasta file.',
                        dest='contig_fasta',
                        default=None)
    parser.add_argument('-d','--input-depth',
                        required=True,
                        nargs='*',
                        help='Path to the input depth file(every position depth generated from mosdepth or bedtools genomecov).If mulptile samples binning ,                        you can input multiple files.',
                        dest='contig_depth',
                        default=None)
    parser.add_argument('-c','--cannot-link',
                        required=True,
                        help='Path to the input can not link file generated from other additional biological information,one row for one can not link                               constraint.The file format:contig_1\tcontig_2.',
                        dest='cannot_link',
                        default=None)
    parser.add_argument('-o','--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None)
    parser.add_argument('--checkm-posprocess',
                        required=False,
                        type=bool,
                        help='Use CheckM to do the postprocess of the bins or not',
                        dest='checkm',
                        default=True
                        )

    return parser.parse_args()

def validate_args(args):
    def except_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
                sys.exit(1)

    def except_file_depth(f):
        if f is not None:
            for depth_path in f:
                if not os.path.exists(depth_path):
                    sys.stderr.write(f"Error: Expected file '{depth_path}' does not exist\n")
                    sys.exit(1)
    except_file(args.contig_fasta)
    except_file_depth(args.contig_depth)
    except_file(args.cannot_link)

def get_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    basepair_sum = 0
    threshold = 0
    whole_len = np.sum(contig_len)
    contig_len.sort(reverse = True)
    index = 0
    while(basepair_sum / whole_len < 0.95):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    threshold = max(threshold , 5000)
    return threshold

def main(args=None):
    if args is None:
        args = sys.argv

    args = parse_args(args)
    validate_args(args)

    logger = logging.getLogger('Semi-MetaBin')
    logger.setLevel(logging.INFO)

    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    out = args.output
    if not os.path.exists(out):
        os.mkdir(out)

    #Processing input contig fasta file
    logger.info('Start generating kmer features from fasta file with threshold length 1000.')

    kmer_1000 = generate_kmer_features_from_fasta(args.contig_fasta,1000,4)

    short_contig_bp = 0
    whole_contig_bp = 0
    contig_length_list = []

    for seq_record in SeqIO.parse(args.contig_fasta , "fasta"):
        contig_length_list.append(len(seq_record))
        whole_contig_bp += len(seq_record)
        if len(seq_record) > 1000 and len(seq_record) <= 1500:
            short_contig_bp += len(seq_record)

    # if considering contigs between 1000bp and 1500 bp in binning
    if short_contig_bp / whole_contig_bp > 0.05:
        is_short_binned = True
    else:
        is_short_binned = False

    # threshold for generating must link pairs
    threshold = get_threshold(contig_length_list)

    if not is_short_binned:
        logger.info('Generating kmer features from fasta file with threshold length 1500(do not bin contigs shorter than 1500).')
        kmer_1500 = generate_kmer_features_from_fasta(args.contig_fasta, 1500, 4)

    logger.info('Generating kmer features for must link pair.')
    kmer_split = generate_kmer_features_from_fasta(args.contig_fasta,1000,4,split=True,threshold=threshold)

    logger.info('Calculating coverage for every sample.')

    # generating coverage for every contig and for must link pair
    logger.info('Processing Sample{}'.format(1))
    contig_cov , must_link_contig_cov = calculate_coverage(args.contig_depth[0],threshold)

    if len(args.contig_depth) > 1:
        for index_sample , depth_sample in enumerate(args.contig_depth):
            if index_sample == 0:
                continue
            logger.info('Processing Sample{}'.format(index_sample + 1))
            contig_cov_ , must_link_contig_cov_ = calculate_coverage(depth_sample,threshold)
            contig_cov = pd.merge(contig_cov,contig_cov_,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)
            must_link_contig_cov = pd.merge(must_link_contig_cov,must_link_contig_cov_,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    data_1000 = pd.merge(kmer_1000,contig_cov,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    data_split = pd.merge(kmer_split,must_link_contig_cov,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    if not is_short_binned:
        data_1500 = pd.merge(kmer_1500,contig_cov,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    logger.info('Estimating number of bins using single-copy marker genes.')


    contig_output = os.path.join(args.output,os.path.split(args.contig_fasta)[1] + '.frag')
    if not os.path.exists(contig_output):
        frag_out_log = open(contig_output + '.out','w')
        subprocess.check_call(
            ['run_FragGeneScan.pl',
             '-genome={}'.format(args.contig_fasta),
             '-out={}'.format(contig_output),
             '-complete=0',
             '-train=complete',
             '-thread=48',
             ],
            stdout = frag_out_log,
            stderr = subprocess.DEVNULL,
        )

    hmm_output = os.path.join(args.output,os.path.split(args.contig_fasta)[1] + '.hmmout')
    if not os.path.exists(hmm_output):
        hmm_out_log = open(hmm_output+'.out','w')
        subprocess.check_call(
            ['hmmsearch',
             '--domtblout',
             hmm_output,
             '--cut_tc',
             '--cpu', str(48),
             os.path.split(__file__)[0] + '/marker.hmm',
             contig_output+'.faa',
            ],
            stdout=hmm_out_log,
            stderr=subprocess.DEVNULL,
        )


    seed_1000_output = os.path.join(args.output,os.path.split(args.contig_fasta)[1] + '_1000.seed')
    seed_1500_output = os.path.join(args.output,os.path.split(args.contig_fasta)[1] + '_1500.seed')

    if not os.path.exists(seed_1000_output):
        getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
        subprocess.check_call(
            [   'perl', getmarker,
                hmm_output,
                args.contig_fasta,
                str(1001),
                seed_1000_output,
            ]
        )

    if not is_short_binned:
        if not os.path.exists(seed_1500_output):
            getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
            subprocess.check_call(
                ['perl', getmarker,
                 hmm_output,
                 args.contig_fasta,
                 str(1501),
                 seed_1500_output,
                 ]
            )

    if is_short_binned:
        seed_1000 = open(seed_1000_output).read().split('\n')
        logger.info('The number of bins:', len(seed_1000))
    else:
        seed_1000 = open(seed_1000_output).read().split('\n')
        seed_1500 = open(seed_1500_output).read().split('\n')
        logger.info('The number of bins:', len(seed_1500))

    







if __name__ == '__main__':
    main(sys.argv)
