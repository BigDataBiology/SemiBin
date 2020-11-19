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
from .semi_supervised_model import loss_function,unsupervised_feature_Dataset,feature_Dataset,Semi_encoding
from torch.utils.data import DataLoader
import torch
from torch.optim import lr_scheduler
from tqdm import tqdm
from sklearn.cluster import KMeans
from sklearn.ensemble import IsolationForest
from atomicwrites import atomic_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
    parser.add_argument('--epoches',
                        required=False,
                        type=int,
                        help='Epoches used in the training process.',
                        dest='epoches',
                        default=20
    )
    parser.add_argument('--batch-size',
                        required=False,
                        type=int,
                        help='Batch size used in the training process.',
                        dest= 'batchsize',
                        default=1024,
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
    contig_length_dict = {}
    contig_dict = {}
    for seq_record in SeqIO.parse(args.contig_fasta , "fasta"):
        contig_length_list.append(len(seq_record))
        whole_contig_bp += len(seq_record)
        if len(seq_record) > 1000 and len(seq_record) <= 1500:
            short_contig_bp += len(seq_record)
        contig_length_dict[str(seq_record.id).strip('')] = len((seq_record.seq))
        contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

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

    contig_cov = contig_cov.apply(lambda x: x+1e-5)
    must_link_contig_cov = must_link_contig_cov.apply(lambda x: x+1e-5)
    if len(args.contig_depth) == 1:
        contig_cov = contig_cov / 100
        must_link_contig_cov = must_link_contig_cov / 100
    else:
        contig_cov = contig_cov.div(contig_cov.sum(axis=0), axis=1)
        contig_cov = contig_cov.div(contig_cov.sum(axis=1), axis=0)

        must_link_contig_cov = must_link_contig_cov.div(must_link_contig_cov.sum(axis=0), axis=1)
        must_link_contig_cov = must_link_contig_cov.div(must_link_contig_cov.sum(axis=1), axis=0)

    data_1000 = pd.merge(kmer_1000,contig_cov,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    data_split = pd.merge(kmer_split,must_link_contig_cov,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    if not is_short_binned:
        data_1500 = pd.merge(kmer_1500,contig_cov,how='inner', on=None,
                    left_index=True, right_index=True, sort=False, copy=True)

    logger.info('Estimating number of bins using single-copy marker genes.')


    contig_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '.frag')
    if not os.path.exists(contig_output + '.faa'):
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

    hmm_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '.hmmout')
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


    seed_1000_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '_1000.seed')
    seed_1500_output = os.path.join(out,os.path.split(args.contig_fasta)[1] + '_1500.seed')

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
        seed_1000 =  [contig for contig in seed_1000 if contig != '' ]
        logger.info('The number of bins:{}'.format(len(seed_1000)))
    else:
        seed_1000 = open(seed_1000_output).read().split('\n')
        seed_1000 =  [contig for contig in seed_1000 if contig != '' ]
        seed_1500 = open(seed_1500_output).read().split('\n')
        seed_1500 =  [contig for contig in seed_1500 if contig != '' ]
        logger.info('The number of bins:{}'.format(len(seed_1500)))

    logger.info('Generate training data:')

    # generate two inputs
    train_input_1 = []
    train_input_2 = []
    train_labels = []

    # can not link
    cannot_link = pd.read_csv(args.cannot_link, sep=',',
                              header=None).values

    namelist = data_1000.index.tolist()
    mapObj = dict(zip(namelist, range(len(namelist))))
    train_data = data_1000.values
    train_data_must_link = data_split.values

    # can not link from contig annotation
    for link in cannot_link:
        train_input_1.append(train_data[mapObj[link[0]]])
        train_input_2.append(train_data[mapObj[link[1]]])
        train_labels.append(0)

    #can not link from bin seed
    for i in range(len(seed_1000)):
        for j in range(i+1,len(seed_1000)):
            train_input_1.append(train_data[mapObj[seed_1000[i]]])
            train_input_2.append(train_data[mapObj[seed_1000[j]]])
            train_labels.append(0)

    # must link from breaking up
    for i in range(0, len(train_data_must_link), 2):
        train_input_1.append(train_data_must_link[i])
        train_input_2.append(train_data_must_link[i + 1])
        train_labels.append(1)

    logger.info('Number of must link pair:{}'.format(train_labels.count(1)))
    logger.info('Number of can not link pair:{}'.format(train_labels.count(0)))

    logger.info('Training model...')


    dataset = feature_Dataset(train_input_1,train_input_2,train_labels)
    train_loader = DataLoader(dataset=dataset, batch_size=args.batchsize, shuffle=True, num_workers=16)

    unlabeled_x = train_data
    unlabeled_train_input1 = []
    unlabeled_train_input2 = []
    for i in range(len(unlabeled_x)):
                    unlabeled_train_input1.append(unlabeled_x[i])
                    unlabeled_train_input2.append(unlabeled_x[i])


    dataset_unlabeled = unsupervised_feature_Dataset(unlabeled_train_input1,unlabeled_train_input2)
    train_loader_unlabeled = DataLoader(dataset=dataset_unlabeled, batch_size=args.batchsize, shuffle=True, num_workers=16)

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu")

    model = Semi_encoding(train_data.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)

    for epoch in tqdm(range(args.epoches)):
        for train_input1, train_input2, train_labels in train_loader:
            model.train()
            train_input1 = train_input1.to(device)
            train_input2 = train_input2.to(device)
            train_labels = train_labels.to(device)
            embedding1, embedding2 = model.forward(train_input1.float(), train_input2.float())
            decoder1, decoder2 = model.decoder(embedding1, embedding2)
            optimizer.zero_grad()
            loss, supervised_loss, unsupervised_loss = loss_function(embedding1.double(), embedding2.double(),
                                                                     train_labels.double(), train_input1.double(),
                                                                     train_input2.double(), decoder1.double(),
                                                                     decoder2.double())
            loss = loss.to(device)
            loss.backward()
            optimizer.step()

        for unlabeled_train_input1 , unlabeled_train_input2 in train_loader_unlabeled:
            model.train()
            unlabeled_train_input1 = unlabeled_train_input1.to(device)
            unlabeled_train_input2 = unlabeled_train_input2.to(device)
            embedding1, embedding2 = model.forward(unlabeled_train_input1.float(), unlabeled_train_input2.float())
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

    with torch.no_grad():
        model.eval()
        if is_short_binned:
            namelist = data_1000.index.tolist()
            row_index = data_1000._stat_axis.values.tolist()
            bin_data = data_1000.values
            init_seed = seed_1000
            num_bin = len(seed_1000)
        else:
            namelist = data_1500.index.tolist()
            row_index = data_1500._stat_axis.values.tolist()
            bin_data = data_1500.values
            init_seed = seed_1500
            num_bin = len(seed_1500)

        seed_index = []
        for temp in init_seed:
            seed_index.append(row_index.index(temp))

        length_weight = [contig_length_dict[name] for name in namelist]
        x = torch.from_numpy(bin_data).to(device)
        embedding = model.embedding(x.float()).detach().cpu().numpy()
        seeds_embedding = embedding[seed_index]
        logger.info('Weighted kmeans clustering.')
        kmeans = KMeans(n_clusters=num_bin, init=seeds_embedding, verbose=1)
        kmeans.fit(embedding, sample_weight=length_weight)

        logger.info('Remove outliers in bins with IsolationForest.')
        contig_labels = kmeans.labels_
        label_set = set(contig_labels)
        for label in label_set:
            index = contig_labels == label
            x_temp = embedding[index]
            clf = IsolationForest(contamination=0.05)
            iforest_label = clf.fit_predict(x_temp)
            iforest_label = np.array([label if temp == 1 else -1 for temp in iforest_label])
            contig_labels[index] = iforest_label


        bin_result = pd.DataFrame(np.concatenate((np.array(namelist).reshape(len(namelist), 1), contig_labels.reshape(len(contig_labels), 1)), axis=1))

        with atomic_write(os.path.join(out, '/bin_result.tsv') , overwrite=True) as ofile:
            bin_result.to_csv(ofile, sep='\t', index=False)

        #write to bins
        res = {}
        for i in range(len(namelist)):
            if contig_labels[i] not in res and contig_labels[i] != -1:
                res[contig_labels[i]] = []
                res[contig_labels[i]].append(namelist[i])
            if contig_labels[i] in res:
                res[contig_labels[i]].append(namelist[i])

        if not os.path.exists(os.path.join(out,'output_bins')):
            os.mkdir(os.path.join(out,'output_bins'))

        for label in res:
            bin = []
            for contig in res[label]:
                rec = SeqRecord(Seq(str(contig_dict[contig])), id=contig, description='')
                bin.append(rec)
            with atomic_write(os.path.join(out,'output_bins','bin.{}.fa'.format(label)), overwrite=True) as ofile:
                SeqIO.write(bin, ofile, 'fasta')

        logger.info('Binning finished.')







if __name__ == '__main__':
    main(sys.argv)
