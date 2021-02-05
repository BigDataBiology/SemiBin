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
    parser.add_argument('-i','--input-fasta',
                        required=True,
                        help='Path to the input contig fasta file.',
                        dest='contig_fasta',
                        default=None)
    parser.add_argument('-d','--input-depth',
                        required=True,
                        nargs='*',
                        help='Path to the input depth file (every position depth generated from mosdepth or bedtools genomecov). '
                            'If using multiple sample binning , you can input multiple files.',
                        dest='contig_depth',
                        default=None)
    parser.add_argument('-c','--cannot-link',
                        required=True,
                        help='Path to the input can not link file generated from other additional biological information, '
                             'one row for each cannot link constraint. '
                             'The file format: contig_1<TAB>contig_2.',
                        dest='cannot_link',
                        default=None)
    parser.add_argument('-o','--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None)
    parser.add_argument('--epoches',
                        required=False,
                        type=int,
                        help='Number of epoches used in the training process.',
                        dest='epoches',
                        default=20
    )
    parser.add_argument('--batch-size',
                        required=False,
                        type=int,
                        help='Batch size used in the training process.',
                        dest='batchsize',
                        default=2048,
                        )
    parser.add_argument('--max-edges',
                        required=False,
                        type=int,
                        help='The maximun number of edges that can be connected to one contig.',
                        dest='max_edges',
                        default=200)
    parser.add_argument('--max-node',
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
    while(basepair_sum / whole_len < 0.98):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    threshold = max(threshold, 4000)
    return threshold

def write_bins(namelist, contig_labels, output, contig_dict, recluster=False, origin_label=0):
    from collections import defaultdict
    res = defaultdict(list)
    for label,name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    for label in res:
        bin = []
        whole_bin_bp = 0
        for contig in res[label]:
            rec = SeqRecord(Seq(str(contig_dict[contig])), id=contig, description='')
            bin.append(rec)
            whole_bin_bp += len(str(contig_dict[contig]))
        if not recluster:
            if whole_bin_bp >= 200000:
                with atomic_write(os.path.join(output, 'bin.{}.fa'.format(label)), overwrite=True) as ofile:
                    SeqIO.write(bin, ofile, 'fasta')
        else:
            if whole_bin_bp >= 200000:
                with atomic_write(os.path.join(output, 'recluster_{0}.bin.{1}.fa'.format(origin_label,label)), overwrite=True) as ofile:
                    SeqIO.write(bin, ofile, 'fasta')


def cal_kl(m1,m2,v1,v2):
        m1 = max(m1,1e-6)
        m2 = max(m2,1e-6)
        v1 = 1 if v1 < 1 else v1
        v2 = 1 if v2 < 1 else v2
        value = np.log(np.sqrt(v2 / v1)) + np.divide(np.add(v1,np.square(m1 - m2)),2 * v2) - 0.5
        return min(max(value,1e-6),1-1e-6)

def cal_num_bins(fasta_path,contig_output,hmm_output,seed_output,binned_short):
    if not os.path.exists(contig_output + '.faa'):
        frag_out_log = open(contig_output + '.out','w')
        subprocess.check_call(
            ['run_FragGeneScan.pl',
             '-genome={}'.format(fasta_path),
             '-out={}'.format(contig_output),
             '-complete=0',
             '-train=complete',
             '-thread=48',
             ],
            stdout = frag_out_log,
            stderr = subprocess.DEVNULL,
        )


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




    if not os.path.exists(seed_output):
            if binned_short:
                getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
                subprocess.check_call(
                    ['perl', getmarker,
                     hmm_output,
                     fasta_path,
                     str(1001),
                     seed_output,
                     ],
                    stderr=subprocess.DEVNULL,
                )
            else:
                getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
                subprocess.check_call(
                    ['perl', getmarker,
                     hmm_output,
                     fasta_path,
                     str(2501),
                     seed_output,
                     ],
                    stderr=subprocess.DEVNULL,
                )




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
    os.makedirs(out, exist_ok=True)

    whole_contig_bp = 0
    contig_bp_2500 = 0
    contig_length_list = []
    contig_length_dict = {}
    contig_dict = {}
    for seq_record in SeqIO.parse(args.contig_fasta , "fasta"):
        if len(seq_record) >= 1000 and len(seq_record) <= 2500:
            contig_bp_2500 += len(seq_record)
        contig_length_list.append(len(seq_record))
        whole_contig_bp += len(seq_record)
        contig_length_dict[str(seq_record.id).strip('')] = len((seq_record.seq))
        contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

    # threshold for generating must link pairs
    threshold = get_threshold(contig_length_list)

    binned_short = contig_bp_2500 / whole_contig_bp < 0.05

    #Processing input contig fasta file
    logger.info('Start generating kmer features from fasta file.')
    #if not os.path.exists(os.path.join(out, 'data.csv')):

    kmer_whole = generate_kmer_features_from_fasta(args.contig_fasta,
            1000 if binned_short else 2500,
            4)

    logger.info('Generating kmer features for must link pair.')
    #if not os.path.exists(os.path.join(out,'data_split.csv')):
    kmer_split = generate_kmer_features_from_fasta(args.contig_fasta,1000,4,split=True,threshold=threshold)

    logger.info('Calculating coverage for every sample.')

    # generating coverage for every contig and for must link pair
    n_sample = len(args.contig_depth)

    is_combined = n_sample >= 5

    if is_combined:
        logger.info('Processing Sample{}'.format(1))
        contig_cov , must_link_contig_cov = calculate_coverage(args.contig_depth[0],threshold,is_combined=is_combined)
        if len(args.contig_depth) > 1:
            for index_sample , depth_sample in enumerate(args.contig_depth):
                if index_sample == 0:
                    continue
                logger.info('Processing Sample{}'.format(index_sample + 1))
                contig_cov_ , must_link_contig_cov_ = calculate_coverage(depth_sample,threshold,is_combined=is_combined)
                contig_cov = pd.merge(contig_cov,contig_cov_,how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)
                must_link_contig_cov = pd.merge(must_link_contig_cov,must_link_contig_cov_,how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)

        contig_cov = contig_cov.apply(lambda x: x+1e-5)
        must_link_contig_cov = must_link_contig_cov.apply(lambda x: x+1e-5)
        contig_cov = contig_cov / 100
        must_link_contig_cov = must_link_contig_cov / 100
        data = pd.merge(kmer_whole, contig_cov, how='inner', on=None,
                             left_index=True, right_index=True, sort=False, copy=True)

        data_split = pd.merge(kmer_split, must_link_contig_cov, how='inner', on=None,
                              left_index=True, right_index=True, sort=False, copy=True)
    else:
        logger.info('Processing Sample{}'.format(1))
        contig_cov = calculate_coverage(args.contig_depth[0],threshold,is_combined=is_combined)
        if len(args.contig_depth) > 1:
            for index_sample , depth_sample in enumerate(args.contig_depth):
                if index_sample == 0:
                    continue
                logger.info('Processing Sample{}'.format(index_sample + 1))
                contig_cov_  = calculate_coverage(depth_sample,threshold,is_combined=is_combined)
                contig_cov = pd.merge(contig_cov,contig_cov_,how='inner', on=None,
                        left_index=True, right_index=True, sort=False, copy=True)

        contig_cov = contig_cov.apply(lambda x: x+1e-5)

        data = pd.merge(kmer_whole, contig_cov, how='inner', on=None,
                             left_index=True, right_index=True, sort=False, copy=True)
        data_split = kmer_split

    # data = pd.merge(kmer_whole, contig_cov, how='inner', on=None,
    #                      left_index=True, right_index=True, sort=False, copy=True)
    # data_split = kmer_split

    with atomic_write(os.path.join(out, 'data.csv'), overwrite=True) as ofile:
        data.to_csv(ofile)

    with atomic_write(os.path.join(out,'data_split.csv'), overwrite=True) as ofile:
        data_split.to_csv(ofile)


    # data = pd.read_csv(os.path.join(out, 'data.csv'), index_col=0)
    # data_split = pd.read_csv(os.path.join(out, 'data_split.csv'), index_col=0)

    contig_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '.frag')
    hmm_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '.hmmout')
    seed_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '.seed')



    cal_num_bins(args.contig_fasta,contig_output,hmm_output,seed_output,binned_short)

    seed = open(seed_output).read().split('\n')
    seed =  [contig for contig in seed if contig != '' ]


    logger.info('Generate training data:')

    kmer = data.values[:,0:136]
    depth = data.values[:,136:len(data.values[0])]


    # #generate two inputs
    train_input_1 = []
    train_input_2 = []
    train_labels = []

    # can not link
    cannot_link = pd.read_csv(args.cannot_link, sep=',',
                              header=None).values



    namelist = data.index.tolist()
    mapObj = dict(zip(namelist, range(len(namelist))))
    row_index = data._stat_axis.values.tolist()
    train_data = data.values
    train_data_must_link = data_split.values

    if not is_combined:
        train_data_input = train_data[:,0:136]
        train_data_split_input = train_data_must_link
    else:
        train_data_input = train_data
        train_data_split_input = train_data_must_link

    # can not link from contig annotation
    for link in cannot_link:
        train_input_1.append(train_data_input[mapObj[link[0]]])
        train_input_2.append(train_data_input[mapObj[link[1]]])
        train_labels.append(0)

    #can not link from bin seed
    for i in range(len(seed)):
        for j in range(i+1,len(seed)):
            train_input_1.append(train_data_input[mapObj[seed[i]]])
            train_input_2.append(train_data_input[mapObj[seed[j]]])
            train_labels.append(0)

    # must link from breaking up
    for i in range(0, len(train_data_must_link), 2):
        train_input_1.append(train_data_split_input[i])
        train_input_2.append(train_data_split_input[i + 1])
        train_labels.append(1)



    logger.info('Number of must link pair:{}'.format(train_labels.count(1)))
    logger.info('Number of can not link pair:{}'.format(train_labels.count(0)))

    logger.info('Training model...')


    dataset = feature_Dataset(train_input_1,train_input_2,train_labels)
    train_loader = DataLoader(dataset=dataset, batch_size=args.batchsize, shuffle=True, num_workers=0)


    unlabeled_x = train_data_input
    unlabeled_train_input1 = []
    unlabeled_train_input2 = []
    for i in range(len(unlabeled_x)):
                    unlabeled_train_input1.append(unlabeled_x[i])
                    unlabeled_train_input2.append(unlabeled_x[i])


    dataset_unlabeled = unsupervised_feature_Dataset(unlabeled_train_input1,unlabeled_train_input2)
    train_loader_unlabeled = DataLoader(dataset=dataset_unlabeled, batch_size=args.batchsize, shuffle=True, num_workers=16)

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu")

    if not is_combined:
        model = Semi_encoding_single(train_data_input.shape[1]).to(device)
    else:
        model = Semi_encoding_multiple(train_data_input.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)

    for epoch in tqdm(range(args.epoches)):
        for train_input1, train_input2, train_label in train_loader:
            model.train()
            train_input1 = train_input1.to(device)
            train_input2 = train_input2.to(device)
            train_label = train_label.to(device)
            embedding1, embedding2 = model.forward(train_input1.float(), train_input2.float())
            decoder1, decoder2 = model.decoder(embedding1, embedding2)
            optimizer.zero_grad()
            loss, supervised_loss, unsupervised_loss = loss_function(embedding1.double(), embedding2.double(),
                                                                     train_label.double(), train_input1.double(),
                                                                     train_input2.double(), decoder1.double(),
                                                                     decoder2.double())
            loss = loss.to(device)
            loss.backward()
            optimizer.step()

        for unlabeled_train_input1 , unlabeled_train_input2 in train_loader_unlabeled:
            model.train()
            unlabeled_train_input1 = unlabeled_train_input1.to(device)
            unlabeled_train_input2 = unlabeled_train_input2.to(device)
            embedding1, embedding2 = model.forward(unlabeled_train_input1.float(), unlabeled_train_input1.float())
            decoder1, decoder2 = model.decoder(embedding1, embedding2)
            optimizer.zero_grad()
            loss = loss_function(embedding1.double(), embedding2.double(),
                                                                     None, unlabeled_train_input1.double(),
                                                                     unlabeled_train_input2.double(), decoder1.double(),
                                                                     decoder2.double(),is_label=False).to(device)
            loss.backward()
            optimizer.step()
        scheduler.step()


    logger.info('Training finished.')
    logger.info('Start binning...')
    torch.save(model,os.path.join(out,'model.h5'))
    #model = torch.load(os.path.join(out,'model.h5'), map_location=torch.device('cpu'))
    #model = torch.load(os.path.join(out, 'model.h5'))
    #model.to(device)
    with torch.no_grad():
        model.eval()
        x = torch.from_numpy(train_data_input).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()
        embedding_matrix = kneighbors_graph(embedding, n_neighbors=args.max_edges, mode='distance', p=2, n_jobs=-1).toarray()
        kmer_matrix = kneighbors_graph(train_data_input,n_neighbors=args.max_edges, mode='distance', p =2, n_jobs=-1).toarray()
        embedding_matrix[kmer_matrix == 0] = 0

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
                            temp_depth  += 1 - cal_kl(depth[i][2*k], depth[j][2*k], depth[i][2*k+1], depth[j][2*k+1])
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
    result = g.community_infomap(edge_weights=edges_weight,vertex_weights=length_weight)
    contig_labels = np.zeros(shape=(len(matrix)), dtype=np.int)

    for i in range(len(result)):
        temp = result[i]
        for infomap_index in temp:
            contig_labels[infomap_index] = i

    output_bin_path = os.path.join(out,'output_bins')
    if not os.path.exists(output_bin_path):
        os.mkdir(output_bin_path)

    write_bins(namelist, contig_labels, output_bin_path, contig_dict)
    if not is_combined:
        mean_index = [2 * temp for temp in range(n_sample)]
        depth_mean = depth[:, mean_index] / 100
        scaling = np.mean(np.abs(embedding)) / np.mean(depth_mean)
        base = 10
        weight = 2 * base * math.ceil(scaling / base)
        embedding_new = np.concatenate((embedding, depth_mean * weight), axis=1)
    else:
        embedding_new = embedding

    bin_files = os.listdir(output_bin_path)
    logger.info('Reclustering.')

    for bin in bin_files:
        if os.path.exists(os.path.join(output_bin_path, bin)):
            contig_list = []
            for seq_record in SeqIO.parse(os.path.join(output_bin_path, bin), "fasta"):
                contig_list.append(seq_record.id)
            contig_output = os.path.join(output_bin_path, bin) + '.frag'
            hmm_output = os.path.join(output_bin_path, bin) + '.hmmout'
            seed_output = os.path.join(output_bin_path, bin) + '.seed'
            try:
                cal_num_bins(os.path.join(output_bin_path, bin),contig_output,hmm_output,seed_output,binned_short)
            except:
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
                length_weight = np.array([contig_length_dict[name] for name in contig_list])
                seeds_embedding = embedding_new[seed_index]
                kmeans = KMeans(n_clusters=num_bin, init=seeds_embedding,n_init=1)
                kmeans.fit(re_bin_features, sample_weight=length_weight)
                labels = kmeans.labels_
                write_bins(contig_list, labels, os.path.join(out, 'output_recluster_bins'), contig_dict,
                           recluster=True, origin_label=int(bin.split('.')[-2]))
            else:
                shutil.copy(os.path.join(output_bin_path, bin), os.path.join(out, 'output_recluster_bins'))

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main(sys.argv)
