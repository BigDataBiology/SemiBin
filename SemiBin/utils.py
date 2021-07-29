import os
import subprocess
from atomicwrites import atomic_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import random
import shutil

def validate_args(args):
    def expect_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(
                    f"Error: Expected file '{f}' does not exist\n")
                sys.exit(1)

    def expect_file_list(fs):
        if fs is not None:
            for f in fs:
                expect_file(f)

    if args.cmd == 'predict_taxonomy':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)

    if args.cmd == 'generate_data_single':
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'generate_data_multi':
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'train':
        if args.mode == 'single':
            if len(args.contig_fasta) > 1:
                sys.stderr.write(
                    f"Error: Expected one fasta file with single mode.\n")
                sys.exit(1)

            if len(args.data) > 1:
                sys.stderr.write(
                    f"Error: Expected one data.csv file with single mode.\n")
                sys.exit(1)

            if len(args.data_split) > 1:
                sys.stderr.write(
                    f"Error: Expected one data_split.csv file with single mode.\n")
                sys.exit(1)

            if len(args.cannot_link) > 1:
                sys.stderr.write(
                    f"Error: Expected one cannot.txt file with single mode.\n")
                sys.exit(1)

            if args.bams is None:
                sys.stderr.write(
                    f"Error: Need to input bams used with single mode.\n")
                sys.exit(1)

            expect_file(args.contig_fasta[0])
            expect_file_list(args.bams)
            expect_file(args.data[0])
            expect_file(args.data_split[0])
            expect_file(args.cannot_link[0])

        elif args.mode == 'several':
            assert len(args.contig_fasta) == len(args.data) == len(args.data_split) == len(args.cannot_link), 'Must input same number of fasta, data, data_split, cannot files!'

            expect_file_list(args.contig_fasta)
            expect_file_list(args.data)
            expect_file_list(args.data_split)
            expect_file_list(args.cannot_link)

        else:
            sys.stderr.write(
                f"Error: Please use training mode with [single/several].\n")
            sys.exit(1)

    if args.cmd == 'bin':
        if args.environment is None and args.model_path is None:
            sys.stderr.write(
                f"Error: Please choose input a model path or use our built-in model for [human_gut/dog_gut/ocean].\n")
            sys.exit(1)
        if args.environment is not None and args.model_path is not None:
            sys.stderr.write(
                f"Error: Please choose input a model path or use our built-in model for [human_gut/dog_gut/ocean].\n")
            sys.exit(1)
        if args.environment is not None:
            if args.environment not in ['human_gut', 'dog_gut', 'ocean']:
                sys.stderr.write(
                    f"Error: Please choose a built-in model in [human_gut/dog_gut/ocean].\n")
                sys.exit(1)
        if args.model_path is not None:
            expect_file(args.model_path)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)
        expect_file(args.data)

    if args.cmd == 'single_easy_bin':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'multi_easy_bin':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

def get_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    import numpy as np
    basepair_sum = 0
    threshold = 0
    whole_len = np.sum(contig_len)
    contig_len.sort(reverse=True)
    index = 0
    while(basepair_sum / whole_len < 0.98):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    return np.clip(threshold, 4000, None)


def parse_mmseqs(mmseqs_result):
    species_result = mmseqs_result.query('rank_name == "species" and score > 0.95').values
    genus_result = mmseqs_result.query('rank_name == "genus" and score > 0.80').values

    cannot_link_species = []
    for i in range(len(species_result)):
        for j in range(i + 1, len(species_result)):
            if species_result[i][2] != species_result[j][2]:
                cannot_link_species.append(
                    (species_result[i][0], species_result[j][0]))

    cannot_link_genus = []
    for i in range(len(genus_result)):
        for j in range(i + 1, len(genus_result)):
            if genus_result[i][2] != genus_result[j][2]:
                cannot_link_genus.append(
                    (genus_result[i][0], genus_result[j][0]))

    cannot_link_mix = []
    for i in range(len(species_result)):
        genus_name = species_result[i][4].split(';')[-2]
        for j in range(len(genus_result)):
            if genus_name != genus_result[j][2]:
                cannot_link_mix.append(
                    (species_result[i][0], genus_result[j][0]))

    return cannot_link_species, cannot_link_genus, cannot_link_mix


def generate_cannot_link(mmseqs_path,namelist,num_threshold,output,sample):
    import itertools
    import numpy as np
    import pandas as pd
    mmseqs_result = pd.read_csv(mmseqs_path, sep='\t', header=None)
    mmseqs_result.columns = ['contig_name', 'taxon_ID', 'rank_name', 'scientific_name', 'temp_1', 'temp_2', 'temp_3',
                             'score', 'lineage']
    mmseqs_result = mmseqs_result[[
        'contig_name', 'rank_name', 'scientific_name', 'score', 'lineage']]
    mmseqs_result['contig_name'] = mmseqs_result['contig_name'].astype(str)
    mmseqs_result = mmseqs_result[mmseqs_result['contig_name'].isin(
        namelist)]
    cannot_link_species, cannot_link_genus, cannot_link_mix = parse_mmseqs(
        mmseqs_result)
    num_whole_data = np.clip(1000 * num_threshold, None, 4_000_000)
    max_num_mix = int(num_whole_data / 8)
    max_num_genus = int(max_num_mix / 2)
    max_num_species = num_whole_data - max_num_mix - max_num_genus
    if len(cannot_link_mix) > max_num_mix:
        cannot_link_mix = random.sample(cannot_link_mix, max_num_mix)

    if len(cannot_link_genus) > max_num_genus:
        cannot_link_genus = random.sample(cannot_link_genus, max_num_genus)

    if len(cannot_link_species) > max_num_species:
        cannot_link_species = random.sample(cannot_link_species, max_num_species)

    with open(os.path.join(output, f'{sample}.txt'), 'w') as out_cannot_link:
        for cannot in itertools.chain(
                cannot_link_species,
                cannot_link_genus,
                cannot_link_mix):
            out_cannot_link.write(f'{cannot[0]},{cannot[1]}\n')

def cal_num_bins(fasta_path, contig_output, hmm_output,
                 seed_output, binned_length, num_process):
    if not os.path.exists(contig_output + '.faa'):
        frag_out_log = open(contig_output + '.out', 'w')
        subprocess.call(
            [shutil.which('FragGeneScan'),
             '-s', fasta_path,
             '-o', contig_output,
             '-w', str(0),
             '-t', 'complete',
             '-p', str(num_process),
             ],
            stdout=frag_out_log,
        )

    if not os.path.exists(hmm_output):
        hmm_out_log = open(hmm_output + '.out', 'w')
        subprocess.check_call(
            ['hmmsearch',
             '--domtblout',
             hmm_output,
             '--cut_tc',
             '--cpu', str(num_process),
             os.path.split(__file__)[0] + '/marker.hmm',
             contig_output + '.faa',
             ],
            stdout=hmm_out_log,
        )

    if not os.path.exists(seed_output):
        getmarker = os.path.split(__file__)[0] + '/getmarker.pl'
        subprocess.check_call(
            ['perl', getmarker,
             hmm_output,
             fasta_path,
             str(binned_length + 1), # threshold
             seed_output,
             ],
            stdout=subprocess.DEVNULL,
        )


def write_bins(namelist, contig_labels, output, contig_dict,
               recluster=False, origin_label=0,minfasta = 200000):
    from collections import defaultdict
    res = defaultdict(list)
    for label, name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    for label in res:
        bin = []
        whole_bin_bp = 0
        for contig in res[label]:
            rec = SeqRecord(
                Seq(str(contig_dict[contig])), id=contig, description='')
            bin.append(rec)
            whole_bin_bp += len(str(contig_dict[contig]))

        ofname = f'bin.{label}.fa' if not recluster \
                    else f'recluster_{origin_label}.bin.{label}.fa'
        if whole_bin_bp >= minfasta:
            with atomic_write(os.path.join(output, ofname), overwrite=True) as ofile:
                SeqIO.write(bin, ofile, 'fasta')


def cal_kl(m, v):
    """
    Calculate KL divergence
    """
    import numpy as np
    m = np.clip(m, 1e-6, None)
    v = np.clip(v, 1.0, None)
    m1 = m.reshape(1, len(m))
    m2 = m.reshape(len(m), 1)

    v1 = v.reshape(1, len(v))
    v2 = v.reshape(len(v), 1)

    value = np.log(np.sqrt(v1 / v2)) +  np.divide(np.square(m1 - m2) + v2, 2 * v1) - 0.5

    return np.clip(value, 1e-6, 1 - 1e-6)

def get_file_md5(fname):
    """
    Calculate Md5 for downloaded file
    """
    import hashlib
    m = hashlib.md5()
    with open(fname,'rb') as fobj:
        while True:
            data = fobj.read(4096)
            if not data:
                break
            m.update(data)

    return m.hexdigest()

def download(logger, GTDB_path):
    """
    Download GTDB.
    GTDB_path: defalt path is $HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB
    """
    import requests
    import tarfile
    logger.info('Downloading GTDB.  It will take a while..')
    GTDB_dir = GTDB_path
    os.makedirs(GTDB_dir, exist_ok=True)

    download_url = 'https://zenodo.org/record/4751564/files/GTDB_v95.tar.gz?download=1'
    download_path = os.path.join(GTDB_dir, 'GTDB_v95.tar.gz')

    with requests.get(download_url, stream=True) as r:
        with open(download_path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    logger.info('Download finished. Checking MD5...')
    if get_file_md5(download_path) == '4a70301c54104e87d5615e3f2383c8b5':
        try:
            tar = tarfile.open(download_path, "r:gz")
            file_names = tar.getnames()
            for file_name in file_names:
                tar.extract(file_name, GTDB_dir)
            tar.close()
        except Exception:
            sys.stderr.write(
                f"Error: cannot unzip the file.")
            sys.exit(1)

        os.remove(download_path)
    else:
        os.remove(download_path)
        sys.stderr.write(
            f"Error: MD5 check failed, removing '{download_path}'.\n")
        sys.exit(1)


def set_random_seed(seed):
    import torch
    import numpy as np
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    np.random.seed(seed)

def process_fasta(fasta_path, ratio):
    """
    Return contig length, contig seq
    """
    whole_contig_bp = 0
    contig_bp_2500 = 0
    contig_length_list = []
    contig_length_dict = {}
    contig_dict = {}

    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        if len(seq_record) >= 1000 and len(seq_record) <= 2500:
            contig_bp_2500 += len(seq_record)
        contig_length_list.append(len(seq_record))
        whole_contig_bp += len(seq_record)
        contig_length_dict[str(seq_record.id).strip(
            '')] = len((seq_record.seq))
        contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

    binned_short = contig_bp_2500 / whole_contig_bp < ratio
    must_link_threshold = get_threshold(contig_length_list)
    return binned_short, must_link_threshold, contig_length_dict, contig_dict

def unzip_fasta(suffix, contig_path):
    import gzip
    import bz2
    if suffix == 'gz':
        contig_name = contig_path.replace(".gz", "")
        ungz_file = gzip.GzipFile(contig_path)
        open(contig_name, "wb+").write(ungz_file.read())
        ungz_file.close()
        return contig_name

    if suffix == 'bz2':
        contig_name = contig_path.replace(".bz2", "")
        unbz2_file = bz2.BZ2File(contig_path)
        open(contig_name, "wb+").write(unbz2_file.read())
        unbz2_file.close()
        return contig_name

def split_data(data, sample, separator):
    """
    split data according their sample in multi-sample binning
    """
    import numpy as np
    part_data = data[data['contig_name'].str.contains(
        '{}'.format(sample + separator))]
    part_data = part_data.set_index('contig_name')
    part_data.index.name = None
    index_list = part_data.index.tolist()
    index_list = [temp.split(separator)[1] for temp in index_list]
    part_data.index = index_list
    abun_scale = (part_data.mean() / 100).apply(np.ceil) * 100
    part_data = part_data.div(abun_scale)

    return part_data

def get_model_path(env):
    if env == 'human_gut':
        model_path = os.path.join(os.path.split(__file__)[0], 'human_gut_model.h5')
        return model_path
    if env == 'dog_gut':
        model_path = os.path.join(os.path.split(__file__)[0], 'dog_gut_model.h5')
        return model_path
    if env == 'ocean':
        model_path = os.path.join(os.path.split(__file__)[0], 'ocean_model.h5')
        return model_path
