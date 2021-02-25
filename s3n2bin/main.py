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
from .semi_supervised_model import loss_function,unsupervised_feature_Dataset,feature_Dataset,Semi_encoding_single,Semi_encoding_multiple
from torch.utils.data import DataLoader
import torch
from torch.optim import lr_scheduler
from tqdm import tqdm
from sklearn.cluster import KMeans
from atomicwrites import atomic_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import shutil
from sklearn.neighbors import kneighbors_graph
from igraph import Graph
import warnings




def parse_args(args):


    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural network for metagenomic binning')

    basic = parser.add_argument_group(title='Basic commands', description=None)
    basic.add_argument('-i','--input-fasta',
                        required=True,
                        help='Path to the input contig fasta file.',
                        dest='contig_fasta',
                        default=None)
    basic.add_argument('-b','--input-bam',
                        required=True,
                        nargs='*',
                        help='Path to the input bam file.'
                             'If using mulptile samples binning, you can input multiple files.',
                        dest='bams',
                        default=None)
    basic.add_argument('-c','--cannot-link',
                        required=True,
                        help='Path to the input can not link file generated from other additional biological information,'
                             'one row for each can not link constraint.'
                             'The file format: contig_1<TAB>contig_2.',
                        dest='cannot_link',
                        default=None)
    basic.add_argument('-o','--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None)
    basic.add_argument('-s','--separator',
                        required=False,
                        type=str,
                        help='Used when multiple samples binning to separete sample name and contig name.(None means single sample and coassemble binning)',
                        dest='separator',
                        default=None,
                       )

    optional = parser.add_argument_group(title='Optional commands', description=None)
    optional.add_argument('-p','--processes',
                        required=False,
                        type=int,
                        help='Number of subprocess used in processing bam files()',
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
                        help='The maximun number of edges that can be connected to one contig.',
                        dest='max_edges',
                        default=200)
    optional.add_argument('--max-node',
                        required=False,
                        type=float,
                        dest='max_node',
                        default=1,
                        help='Percentage of contigs that considered to be binned.')

    return parser.parse_args()

def validate_args(args):

    def except_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
                sys.exit(1)

    def except_file_bam(f):
        if f is not None:
            for bam_path in f:
                if not os.path.exists(bam_path):
                    sys.stderr.write(f"Error: Expected file '{bam_path}' does not exist\n")
                    sys.exit(1)

    except_file(args.contig_fasta)
    except_file_bam(args.bams)
    except_file(args.cannot_link)


def main(args=None):
    if args is None:
        args = sys.argv

    args = parse_args(args)
    validate_args(args)


    logger = logging.getLogger('S3N2Bin')
    logger.setLevel(logging.INFO)

    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    out = args.output
    os.makedirs(out, exist_ok=True)

    multi_binning = False
    if args.separator is not None:
        multi_binning = True



if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main(sys.argv)
