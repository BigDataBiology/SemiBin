import os
import subprocess
import contextlib
import multiprocessing
from atomicwrites import atomic_write
import tempfile
import sys
import random
import shutil

from .fasta import fasta_iter

def validate_normalize_args(logger, args):
    '''Validate and normalize command line arguments'''
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

    if args.cmd != 'download_GTDB' and args.cmd != 'check_install' and args.cmd != 'concatenate_fasta':
        if args.num_process == 0:
            args.num_process = multiprocessing.cpu_count()
            logger.info(f'Setting number of CPUs to {args.num_process}')

        if args.num_process > multiprocessing.cpu_count():
            args.num_process = multiprocessing.cpu_count()

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train', 'bin']:
        if args.orf_finder not in ['prodigal', 'fraggenescan']:
            sys.stderr.write(
                f"Error: SemiBin only supports prodigal or fraggenescan as the ORF finder (--orf_finder option).\n")
            sys.exit(1)
        if args.engine not in ['auto', 'gpu', 'cpu']:
            sys.stderr.write(
                f"Error: You need to specify the engine in auto/gpu/cpu.\n")
            sys.exit(1)

    if args.cmd == 'generate_cannot_links':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.taxonomy_results_fname)

    if args.cmd == 'generate_sequence_features_single':
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'generate_sequence_features_multi':
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
            if args.environment not in ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'global']:
                sys.stderr.write(
                    f"Error: Please choose a built-in model in [human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global].\n")
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

        if args.environment is not None:
            if args.environment not in ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'global']:
                sys.stderr.write(
                    f"Error: Please choose a built-in model in [human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global].\n")
                sys.exit(1)

    if args.cmd == 'multi_easy_bin':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'concatenate_fasta':
        expect_file_list(args.contig_fasta)
        contig_name = []
        for contig in args.contig_fasta:
            contig_name.append(os.path.basename(contig).split('.')[0])
        if len(set(contig_name)) != len(contig_name):
            sys.stderr.write(
                f"Error: Make sure that every contig file have different names.\n")
            sys.exit(1)



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


def generate_cannot_link(mmseqs_path, namelist, num_threshold, output,sample):
    import itertools
    import numpy as np
    import pandas as pd
    mmseqs_result = pd.read_csv(mmseqs_path,
                        sep='\t',
                        header=None,
                        names=['contig_name', 'taxon_ID', 'rank_name',
                                'scientific_name', 'temp_1', 'temp_2', 'temp_3', 'score', 'lineage'],
                        usecols=['contig_name', 'rank_name', 'scientific_name', 'score', 'lineage'])

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

def get_marker(hmmout, fasta_path=None, min_contig_len=None, multi_mode=False, orf_finder = None):
    import pandas as pd
    data = pd.read_table(hmmout, sep=r'\s+',  comment='#', header=None,
                         usecols=(0,3,5,15,16), names=['orf', 'gene', 'qlen', 'qstart', 'qend'])
    if not len(data):
        return []
    data['gene'] = data['gene'].map(lambda m: normalize_marker_trans__dict.get(m , m))
    qlen = data[['gene','qlen']].drop_duplicates().set_index('gene')['qlen']

    def contig_name(ell):
        if orf_finder == 'prodigal':
            contig,_ = ell.rsplit( '_', 1)
        else:
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

def prodigal(contig_file, contig_output):
        with open(contig_output + '.out', 'w') as prodigal_out_log:
            subprocess.check_call(
                ['prodigal',
                 '-i', contig_file,
                 '-p', 'meta',
                 '-q',
                 '-m', # See https://github.com/BigDataBiology/SemiBin/issues/87
                 '-a', contig_output
                 ],
                stdout=prodigal_out_log,
            )

def run_prodigal(fasta_path, num_process, output):
    from .error import LoggingPool

    contigs = {}
    for h, seq in fasta_iter(fasta_path):
        contigs[h] = seq

    total_len = sum(len(s) for s in contigs.values())
    split_len = total_len // num_process

    cur = split_len + 1
    next_ix = 0
    out = None
    with contextlib.ExitStack() as stack:
        for h,seq in contigs.items():
            if cur > split_len and next_ix < num_process:
                if out is not None:
                    out.close()
                out = open(os.path.join(output, 'contig_{}.fa'.format(next_ix)), 'wt')
                out = stack.enter_context(out)

                cur = 0
                next_ix += 1
            out.write(f'>{h}\n{seq}\n')
            cur += len(seq)

    with LoggingPool(num_process) if num_process != 0 else LoggingPool() as pool:
        try:
            for index in range(next_ix):
                pool.apply_async(
                    prodigal,
                    args=(
                        os.path.join(output, 'contig_{}.fa'.format(index)),
                        os.path.join(output, 'contig_{}.faa'.format(index)),
                    ))
            pool.close()
            pool.join()
        except:
            sys.stderr.write(
                f"Error: Running prodigal fail\n")
            sys.exit(1)

    contig_output = os.path.join(output, 'contigs.faa')
    with open(contig_output, 'w') as f:
        for index in range(next_ix):
            f.write(open(os.path.join(output, 'contig_{}.faa'.format(index)), 'r').read())
    return contig_output

def run_fraggenescan(fasta_path, num_process, output):
    try:
        contig_output = os.path.join(output, 'contigs.faa')
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
    except:
        sys.stderr.write(
            f"Error: Running fraggenescan failed\n")
        sys.exit(1)
    return contig_output + '.faa'

def cal_num_bins(fasta_path, binned_length, num_process, multi_mode=False, output = None, orf_finder = 'prodigal'):
    '''Estimate number of bins from a FASTA file

    Parameters
    fasta_path: path
    binned_length: int (minimal contig length)
    num_process: int (number of CPUs to use)
    multi_mode: bool, optional (if True, treat input as resulting from concatenating multiple files)
    '''
    with tempfile.TemporaryDirectory() as tdir:
        if output is not None:
            if os.path.exists(os.path.join(output, 'markers.hmmout')):
                return get_marker(os.path.join(output, 'markers.hmmout'), fasta_path, binned_length, multi_mode, orf_finder=orf_finder)
            else:
                os.makedirs(output, exist_ok=True)
                target_dir = output
        else:
            target_dir = tdir

        run_orffinder = run_prodigal if orf_finder == 'prodigal' else run_fraggenescan
        contig_output = run_orffinder(fasta_path, num_process, tdir)

        hmm_output = os.path.join(target_dir, 'markers.hmmout')
        try:
            with open(os.path.join(tdir, 'markers.hmmout.out'), 'w') as hmm_out_log:
                subprocess.check_call(
                    ['hmmsearch',
                     '--domtblout',
                     hmm_output,
                     '--cut_tc',
                     '--cpu', str(num_process),
                     os.path.split(__file__)[0] + '/marker.hmm',
                     contig_output,
                     ],
                    stdout=hmm_out_log,
                )
        except:
            if os.path.exists(hmm_output):
                os.remove(hmm_output)
            sys.stderr.write(
                f"Error: Running hmmsearch fail\n")
            sys.exit(1)

        marker = get_marker(hmm_output, fasta_path, binned_length, multi_mode, orf_finder=orf_finder)

        return marker


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
        sample + separator, regex=False)]
    part_data = part_data.set_index('contig_name')
    part_data.index.name = None
    part_data.index = [ix.split(separator)[1] for ix in part_data.index]
    if is_combined:
        abun_scale = (part_data.mean() / 100).apply(np.ceil) * 100
        part_data = part_data.div(abun_scale)
    return part_data

def get_model_path(env):
    if env in [
            'human_gut',
            'dog_gut',
            'ocean',
            'soil',
            'cat_gut',
            'human_oral',
            'mouse_gut',
            'pig_gut',
            'built_environment',
            'wastewater',
            'global',
            ]:
        return os.path.join(os.path.split(__file__)[0], f'{env}_model.h5')
    else:
        sys.stderr.write(
            f"Error: Expected environment '{env}' does not exist\n")
        sys.exit(1)

def concatenate_fasta(fasta_files, min_length, output, separator):
    cfasta = os.path.join(output, 'concatenated.fa')
    with open(cfasta, 'wt') as concat_out:
        for fasta in fasta_files:
            sample_name = os.path.basename(fasta).split('.')[0]
            for h, seq in fasta_iter(fasta):
                if separator in h:
                    sys.stderr.write(
                        f"Error: the header of the contig '{h}' contains the separator ('{separator}'), please choose another separator.\n")
                    sys.exit(1)
                if len(seq) >= min_length:
                    header = f'{sample_name}{separator}{h}'
                    concat_out.write(f'>{header}\n{seq}\n')

