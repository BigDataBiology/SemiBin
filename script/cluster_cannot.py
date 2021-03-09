import argparse
import os
import sys
import logging
from Bio import SeqIO
from s3n2bin.utils import validate_args, get_threshold
import pandas as pd
from s3n2bin.semi_supervised_model import train
import torch
from atomicwrites import atomic_write
from s3n2bin.cluster import cluster
import warnings
import multiprocessing
from Bio.SeqRecord import SeqRecord
import gzip
import bz2
import shutil

from sklearn.neighbors import kneighbors_graph
from igraph import Graph
import numpy as np
from s3n2bin.utils import cal_kl, write_bins, cal_num_bins
import os
import math
from Bio import SeqIO
from sklearn.cluster import KMeans
import shutil

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural network for metagenomic binning')

    basic = parser.add_argument_group(title='Basic commands', description=None)
    basic.add_argument('-i', '--input-fasta',
                       required=True,
                       help='Path to the input contig fasta file.',
                       dest='contig_fasta',
                       default=None)
    basic.add_argument('-b', '--input-bam',
                       required=True,
                       nargs='*',
                       help='Path to the input bam file.'
                             'If using multiple samples binning, you can input multiple files.',
                       dest='bams',
                       default=None)
    basic.add_argument('-c', '--cannot-link',
                       required=True,
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

    optional = parser.add_argument_group(
        title='Optional commands', description=None)
    optional.add_argument('-p', '--processes',
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

    contig_fasta = args.contig_fasta

    multi_binning = False
    if args.separator is not None:
        multi_binning = True
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

        # threshold for generating must link pairs
        threshold = get_threshold(contig_length_list)
        binned_short = contig_bp_2500 / whole_contig_bp < 0.05
        logger.info('Calculating coverage for every sample.')

        # generating coverage for every contig and for must link pair
        n_sample = len(args.bams)
        is_combined = n_sample >= 5

        bam_list = args.bams


        data = pd.read_csv(os.path.join(out,'data.csv'),index_col=0)
        data_split = pd.read_csv(os.path.join(out,'data_split.csv'),index_col=0)

        train_data = data.values
        if not is_combined:
            train_data_input = train_data[:, 0:136]
        else:
            train_data_input = train_data

        depth = data.values[:, 136:len(data.values[0])]
        namelist = data.index.tolist()
        row_index = data._stat_axis.values.tolist()
        mapObj = dict(zip(namelist, range(len(namelist))))

        x = train_data_input
        embedding = x
        embedding_matrix = kneighbors_graph(
            embedding,
            n_neighbors=args.max_edges,
            mode='distance',
            p=2,
            n_jobs=-1).toarray()
        cannot_link = pd.read_csv(args.cannot_link[0], sep=',', header=None).values
        must_link = pd.read_csv('/share/inspurStorage/home1/pansj/binning/CAMI_high/mmseqs_solidbin/high_must.txt',sep=',',header=None).values

        for temp in cannot_link:
            embedding_matrix[mapObj[temp[0]]][mapObj[temp[1]]] = 2
            embedding_matrix[mapObj[temp[1]]][mapObj[temp[0]]] = 2
        for temp in must_link:
            embedding_matrix[mapObj[temp[0]]][mapObj[temp[1]]] = 1e-6
            embedding_matrix[mapObj[temp[1]]][mapObj[temp[0]]] = 1e-6
        print('must')
        embedding_matrix[embedding_matrix >= 1] = 1
        embedding_matrix[embedding_matrix == 0] = 1
        embedding_matrix = 1 - embedding_matrix

        threshold = 0.95

        while (threshold >= 0):
            temp_threshold = embedding_matrix.copy()
            temp_threshold[temp_threshold <= threshold] = 0
            num = len(list(set(np.where(temp_threshold > 0)[0])))
            if round(num / len(embedding_matrix), 2) >= args.max_node:
                break
            else:
                threshold -= 0.05

        embedding_matrix[embedding_matrix <= threshold] = 0
        if not is_combined:
            logger.info('Calculating depth matrix.')
            depth_matrix = np.zeros(shape=embedding_matrix.shape)
            for i in range(len(embedding_matrix)):
                for j in range(i + 1, len(embedding_matrix)):
                    if embedding_matrix[i][j] > 0:
                        temp_depth = 0
                        for k in range(n_sample):
                            temp_depth += 1 - \
                                          cal_kl(depth[i][2 * k], depth[j][2 * k],
                                                 depth[i][2 * k + 1], depth[j][2 * k + 1])
                        depth_matrix[i][j] = temp_depth / n_sample

            matrix = embedding_matrix * depth_matrix
        else:
            matrix = embedding_matrix

        edges = []
        edges_weight = []

        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):
                if matrix[i][j] > 1e-6:
                    edges.append((i, j))
                    edges_weight.append(matrix[i][j])

        logger.info('Edges:{}'.format(len(edges)))

        g = Graph()
        vertex = list(range(len(matrix)))
        g.add_vertices(vertex)
        g.add_edges(edges)
        length_weight = np.array([contig_length_dict[name] for name in namelist])
        result = g.community_infomap(
            edge_weights=edges_weight,
            vertex_weights=length_weight)
        contig_labels = np.zeros(shape=(len(matrix)), dtype=np.int)

        for i, r in enumerate(result):
            for infomap_index in r:
                contig_labels[infomap_index] = i

        output_bin_path = os.path.join(out, 'output_bins')
        os.makedirs(output_bin_path, exist_ok=True)

        write_bins(namelist, contig_labels, output_bin_path, contig_dict)
        if not is_combined:
            mean_index = [2 * temp for temp in range(n_sample)]
            depth_mean = depth[:, mean_index] / 100
            scaling = np.mean(np.abs(embedding)) / np.mean(depth_mean)
            base = 10
            weight = 2 * base * math.ceil(scaling / base)
            embedding_new = np.concatenate(
                (embedding, depth_mean * weight), axis=1)
        else:
            embedding_new = embedding

        bin_files = os.listdir(output_bin_path)
        logger.info('Reclustering.')

        for bin in bin_files:
            if os.path.exists(os.path.join(output_bin_path, bin)):
                contig_list = []
                for seq_record in SeqIO.parse(
                        os.path.join(output_bin_path, bin), "fasta"):
                    contig_list.append(seq_record.id)
                contig_output = os.path.join(output_bin_path, bin) + '.frag'
                hmm_output = os.path.join(output_bin_path, bin) + '.hmmout'
                seed_output = os.path.join(output_bin_path, bin) + '.seed'
                try:
                    cal_num_bins(
                        os.path.join(
                            output_bin_path,
                            bin),
                        contig_output,
                        hmm_output,
                        seed_output,
                        binned_short)
                except BaseException:
                    pass
                contig_index = [mapObj[temp] for temp in contig_list]
                re_bin_features = embedding_new[contig_index]
                if not os.path.exists(os.path.join(out, 'output_recluster_bins')):
                    os.mkdir(os.path.join(out, 'output_recluster_bins'))

                if os.path.exists(seed_output):
                    seed = open(seed_output).read().split('\n')
                    seed = [contig for contig in seed if contig != '']
                    init_seed = seed
                    num_bin = len(seed)
                    seed_index = []
                    for temp in init_seed:
                        seed_index.append(row_index.index(temp))
                    length_weight = np.array(
                        [contig_length_dict[name] for name in contig_list])
                    seeds_embedding = embedding_new[seed_index]
                    kmeans = KMeans(
                        n_clusters=num_bin,
                        init=seeds_embedding,
                        n_init=1)
                    kmeans.fit(re_bin_features, sample_weight=length_weight)
                    labels = kmeans.labels_
                    write_bins(contig_list, labels, os.path.join(out, 'output_recluster_bins'), contig_dict,
                               recluster=True, origin_label=int(bin.split('.')[-2]))
                else:
                    shutil.copy(os.path.join(
                        output_bin_path, bin), os.path.join(out, 'output_recluster_bins'))





if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main(sys.argv)



