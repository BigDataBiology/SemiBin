import os
import subprocess
from atomicwrites import atomic_write
import tempfile
import sys
import random
import shutil

from .fasta import fasta_iter

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

            expect_file(args.contig_fasta[0])
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


def get_must_link_threshold(contig_lens):
    """
    calculate the threshold length for must link breaking up
    """
    import numpy as np
    contig_lens = np.array(contig_lens)
    contig_lens.sort()
    frac = np.cumsum(contig_lens)/np.sum(contig_lens)
    ix = np.argmax(frac > 0.02) # argmax finds first True element
    threshold = contig_lens[ix]
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


normalize_marker_trans__dict = {
    'TIGR00388': 'TIGR00389',
    'TIGR00471': 'TIGR00472',
    'TIGR00408': 'TIGR00409',
    'TIGR02386': 'TIGR02387',
}

def get_marker(hmmout, fasta_path=None, min_contig_len=None, multi_mode=False):
    import pandas as pd
    data = pd.read_table(hmmout, sep=r'\s+',  comment='#', header=None,
                         usecols=(0,3,5,15,16), names=['orf', 'gene', 'qlen', 'qstart', 'qend'])
    if not len(data):
        return []
    data['gene'] = data['gene'].map(lambda m: normalize_marker_trans__dict.get(m , m))
    qlen = data[['gene','qlen']].drop_duplicates().set_index('gene')['qlen']

    def contig_name(ell):
        contig,_,_,_ = ell.rsplit( '_', 3)
        return contig

    data = data.query('(qend - qstart) / qlen > 0.4').copy()
    data['contig'] = data['orf'].map(contig_name)
    if min_contig_len is not None:
        contig_len = {h:len(seq) for h,seq in fasta_iter(fasta_path)}
        data = data[data['contig'].map(lambda c: contig_len[c] >= min_contig_len)]
    data = data.drop_duplicates(['gene', 'contig'])

    def extract_seeds(vs, sel):
        vs = vs.sort_values()
        median = vs[len(vs) //2]

        # the original version broke ties by picking the shortest query, so we
        # replicate that here:
        candidates = vs.index[vs == median]
        c = qlen.loc[candidates].idxmin()
        r = list(sel.query('gene == @c')['contig'])
        r.sort()
        return r


    if multi_mode:
        data['bin'] = data['orf'].str.split('.',0, expand=True)[0]
        counts = data.groupby(['bin', 'gene'])['orf'].count()
        res = {}
        for b,vs in counts.groupby(level=0):
            cs = extract_seeds(vs.droplevel(0), data.query('bin == @b', local_dict={'b':b}))
            res[b] = [c.split('.',1)[1] for c in cs]
        return res
    else:
        counts = data.groupby('gene')['orf'].count()
        return extract_seeds(counts, data)

def cal_num_bins(fasta_path, binned_length, num_process, multi_mode=False):
    '''Estimate number of bins from a FASTA file

    Parameters
    fasta_path: path
    binned_length: int (minimal contig length)
    num_process: int (number of CPUs to use)
    multi_mode: bool, optional (if True, treat input as resulting from concatenating multiple files)
    '''
    with tempfile.TemporaryDirectory() as tdir:
        contig_output = os.path.join(tdir, 'contigs.faa')
        with open(contig_output + '.out', 'w') as frag_out_log:
            # We need to call FragGeneScan instead of the Perl wrapper because the
            # Perl wrapper does not handle filepaths correctly if they contain spaces
            # This binary does not handle return codes correctly, though, so we
            # cannot use `check_call`:
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

        hmm_output = os.path.join(tdir, 'markers.hmmout')
        with open(hmm_output + '.out', 'w') as hmm_out_log:
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

        return get_marker(hmm_output, fasta_path, binned_length, multi_mode)


def write_bins(namelist, contig_labels, output, contig_dict,
               recluster=False, origin_label=0, minfasta = 200000):
    '''
    Write binned FASTA files

    Returns: list of files written
    '''
    from collections import defaultdict
    res = defaultdict(list)
    for label, name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    written = []
    for label, contigs in res.items():
        whole_bin_bp = sum(len(contig_dict[contig]) for contig in contigs)

        if whole_bin_bp >= minfasta:
            ofname = f'bin.{label}.fa' if not recluster \
                    else f'recluster_{origin_label}.bin.{label}.fa'
            ofname = os.path.join(output, ofname)
            written.append(ofname)
            with atomic_write(ofname, overwrite=True) as ofile:
                for contig in contigs:
                    ofile.write(f'>{contig}\n{contig_dict[contig]}\n')
    return written



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
    Returns

    binned_short: whether to include short contigs
    must_link_threshold: threshold to break up contigs
    contigs: dictionary ID -> contig sequence
    """
    whole_contig_bp = 0
    contig_bp_2500 = 0
    contig_length_list = []
    contig_dict = {}

    for h, seq in fasta_iter(fasta_path):
        if 1000 <= len(seq) <= 2500:
            contig_bp_2500 += len(seq)
        contig_length_list.append(len(seq))
        whole_contig_bp += len(seq)
        contig_dict[h] = seq

    binned_short = contig_bp_2500 / whole_contig_bp < ratio
    must_link_threshold = get_must_link_threshold(contig_length_list)
    return binned_short, must_link_threshold, contig_dict


def split_data(data, sample, separator, is_combined = True):
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
    if is_combined:
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
