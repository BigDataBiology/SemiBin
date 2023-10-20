import os
import subprocess
import multiprocessing
import tempfile
import sys
import random
import contextlib

from .fasta import fasta_iter

@contextlib.contextmanager
def possibly_compressed_write(filename):
    from .atomicwrite import atomic_write
    if filename.endswith('.gz'):
        import gzip
        mode = 'wb'
        transf = lambda f : gzip.open(f, mode='wt')
    elif filename.endswith('.bz2'):
        import bz2
        mode = 'wb'
        transf = lambda f : bz2.open(f, mode='wt')
    elif filename.endswith('.xz'):
        import lzma
        mode = 'wb'
        transf = lambda f : lzma.open(f, mode='wt')
    else:
        mode = 'wt'
        transf = lambda f : f
    with atomic_write(filename, mode=mode, overwrite=True) as f:
        g = transf(f)
        yield g
        if g is not f:
            g.close()

def check_training_type(logger, args):
    if args.is_semibin2:
        args.training_type = 'auto'

    if args.training_type == 'semi' and args.self_supervised:
        logger.error('Both --training-type=semi and --self-supervised arguments used')
        sys.exit(1)

    if args.training_type == 'self' and args.semi_supervised:
        logger.warning('Both --training-type=self and --semi-supervised arguments used')
        sys.exit(1)

    if not args.self_supervised and not args.semi_supervised:
        if args.training_type == 'self' or args.is_semibin2:
            logger.info(
                f"SemiBin will run in self supervised mode")
            args.training_type = 'self'
        else:
            logger.info(
                f"SemiBin will run in semi supervised mode")
            args.training_type = 'semi'

    elif args.self_supervised and args.semi_supervised:
        if args.is_semibin2:
            logger.warning(
                f'You chose both semi-supervised and self-supervised learning! SemiBin will use self-supervised learning')
            args.training_type = 'self'
        else:
            logger.warning(
                f'You chose both semi-supervised and self-supervised learning! SemiBin will use semi-supervised learning (this may change in the future)')
            args.training_type = 'semi'

    elif args.self_supervised:
        args.training_type = 'self'
    else:
        args.training_type = 'semi'


def validate_normalize_args(logger, args):
    '''Validate and normalize command line arguments'''
    exit_with_error = False

    def expect_file(f):
        nonlocal exit_with_error
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(
                    f"Error: Expected file '{f}' does not exist\n")
                exit_with_error = True

    def expect_file_list(fs):
        if fs is not None:
            for f in fs:
                expect_file(f)

    if hasattr(args, 'num_process'):
        if args.num_process == 0:
            args.num_process = multiprocessing.cpu_count()
            logger.info(f'Setting number of CPUs to {args.num_process}')

        if args.num_process > multiprocessing.cpu_count():
            args.num_process = multiprocessing.cpu_count()
        os.environ['NUMEXPR_NUM_THREADS'] = str(args.num_process)
        os.environ['NUMEXPR_MAX_THREADS'] = str(args.num_process)
        os.environ['OMP_NUM_THREADS'] = str(args.num_process)


    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train', 'bin']:
        if args.orf_finder not in ['prodigal', 'fraggenescan', 'fast-naive']:
            sys.stderr.write(
                f"Error: SemiBin only supports 'prodigal'/'fraggenescan'/'fast-naive' as the ORF finder (--orf_finder option).\n")
            exit_with_error = True
        expect_file(args.prodigal_output_faa)
        if args.prodigal_output_faa is not None:
            args.orf_finder = 'prodigal'

    if hasattr(args, 'engine'):
        args.engine = args.engine.lower()
        if args.engine not in ['auto', 'gpu', 'cpu']:
            sys.stderr.write(
                f"Error: Argument '--engine' needs to be one of auto[default]/gpu/cpu.\n")
            exit_with_error = True

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

    if args.cmd in ['train', 'train_semi', 'train_self']:
        if not args.train_from_many:
            if len(args.data) > 1:
                sys.stderr.write(
                    f"Error: Expected one data.csv file with single mode.\n")
                exit_with_error = True

            if len(args.data_split) > 1:
                sys.stderr.write(
                    f"Error: Expected one data_split.csv file with single mode.\n")
                exit_with_error = True
            if args.cmd in ['train_semi', 'train']:
                if len(args.contig_fasta) > 1:
                    sys.stderr.write(
                        f"Error: Expected one fasta file with single mode.\n")
                    exit_with_error = True

                if len(args.cannot_link) > 1:
                    sys.stderr.write(
                        f"Error: Expected one cannot.txt file with single mode.\n")
                    exit_with_error = True
                expect_file(args.cannot_link[0])
                expect_file(args.contig_fasta[0])

            expect_file(args.data[0])
            expect_file(args.data_split[0])

        else:
            if args.cmd in ['train_semi', 'train']:
                assert len(args.contig_fasta) == len(args.data) == len(args.data_split) == len(args.cannot_link), 'Must input same number of fasta, data, data_split, cannot files!'
                expect_file_list(args.cannot_link)
                expect_file_list(args.contig_fasta)
            else:
                assert len(args.data) == len(args.data_split) , 'Must input same number of data, data_split!'

            expect_file_list(args.data)
            expect_file_list(args.data_split)

    if args.cmd in ['single_easy_bin', 'multi_easy_bin']:
        if args.sequencing_type.lower() in ['short', 'short-read', 'short_reads', 'short-reads', 'short_read']:
            args.sequencing_type = 'short_read'
        elif args.sequencing_type.lower() in ['long', 'long-read', 'long_reads', 'long-reads', 'long_read']:
            args.sequencing_type = 'long_read'
        else:
            sys.stderr.write(
                f"Error: Did not understand sequencing_type argument '{args.sequencing_type}' (should be short_reads or long_reads).\n")
            exit_with_error = True
        logger.info(f'Binning for {args.sequencing_type}')

    if args.cmd == 'bin':
        if args.environment is None and args.model_path is None:
            sys.stderr.write(
                f"Error: Please choose input a model path or use our built-in model for [human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global].\n")
            exit_with_error = True
        if args.environment is not None and args.model_path is not None:
            sys.stderr.write(
                f"Error: Please choose input a model path or use our built-in model for [human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global].\n")
            exit_with_error = True
        if args.environment is not None:
            # This triggers checking that the environment is valid
            get_model_path(args.environment)
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
            # This triggers checking that the environment is valid
            get_model_path(args.environment)
        else:
            check_training_type(logger, args)

    if args.cmd == 'multi_easy_bin':
        check_training_type(logger, args)
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

        if args.training_type not in ['semi', 'self']:
            sys.stderr.write(
                f"Error: You need to specify the training algorithm in semi/self.\n")
            exit_with_error = True

    if args.cmd == 'concatenate_fasta':
        expect_file_list(args.contig_fasta)
        contig_name = []
        for contig in args.contig_fasta:
            contig_name.append(os.path.basename(contig).split('.')[0])
        if len(set(contig_name)) != len(contig_name):
            sys.stderr.write(
                f"Error: Make sure that every contig file have different names.\n")
            exit_with_error = True

    if getattr(args, 'train_from_many', False):
        args.mode = 'several'
    elif args.cmd in ['train', 'train_semi', 'train_self'] and not hasattr(args, 'mode'):
        args.mode = 'single'

    if exit_with_error:
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

def get_marker(hmmout, fasta_path=None, min_contig_len=None, multi_mode=False, orf_finder = None, contig_to_marker = False):
    '''Parse HMM output file and return markers
    '''
    import pandas as pd
    data = pd.read_table(hmmout, sep=r'\s+',  comment='#', header=None,
                         usecols=(0,3,5,15,16), names=['orf', 'gene', 'qlen', 'qstart', 'qend'])
    if not len(data):
        return []
    data['gene'] = data['gene'].map(lambda m: normalize_marker_trans__dict.get(m , m))
    qlen = data[['gene','qlen']].drop_duplicates().set_index('gene')['qlen']

    def contig_name(ell):
        if orf_finder in ['prodigal', 'fast-naive']:
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

    if contig_to_marker:
        from collections import defaultdict
        marker = data['gene'].values
        contig = data['contig'].values
        sequence2markers = defaultdict(list)
        for m, c in zip(marker, contig):
            sequence2markers[c].append(m)
        return sequence2markers
    else:
        def extract_seeds(vs, sel):
            vs = vs.sort_values()
            median = vs.iloc[len(vs) //2]

            # the original version broke ties by picking the shortest query, so we
            # replicate that here:
            candidates = vs.index[vs == median]
            c = qlen.loc[candidates].idxmin()
            r = list(sel.query('gene == @c')['contig'])
            r.sort()
            return r


        if multi_mode:
            data['bin'] = data['orf'].str.split(pat='.', n=0, expand=True)[0]
            counts = data.groupby(['bin', 'gene'])['orf'].count()
            res = {}
            for b,vs in counts.groupby(level=0):
                cs = extract_seeds(vs.droplevel(0), data.query('bin == @b', local_dict={'b':b}))
                res[b] = [c.split('.',1)[1] for c in cs]
            return res
        else:
            counts = data.groupby('gene')['orf'].count()
            return extract_seeds(counts, data)


def cal_num_bins(fasta_path, binned_length, num_process, multi_mode=False, output = None, orf_finder = 'prodigal', prodigal_output_faa=None):
    '''Estimate number of bins from a FASTA file

    Parameters
    fasta_path: path
    binned_length: int (minimal contig length)
    num_process: int (number of CPUs to use)
    multi_mode: bool, optional (if True, treat input as resulting from concatenating multiple files)
    '''
    from .orffinding import run_orffinder
    with tempfile.TemporaryDirectory() as tdir:
        if output is not None:
            if os.path.exists(os.path.join(output, 'markers.hmmout')):
                return get_marker(os.path.join(output, 'markers.hmmout'), fasta_path, binned_length, multi_mode, orf_finder=orf_finder)
            else:
                os.makedirs(output, exist_ok=True)
                target_dir = output
        else:
            target_dir = tdir

        contig_output = run_orffinder(fasta_path, num_process, tdir, orf_finder, prodigal_output_faa=prodigal_output_faa)

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

        return get_marker(hmm_output, fasta_path, binned_length, multi_mode, orf_finder=orf_finder)


def write_bins(namelist, contig_labels, output, contig_seqs,
               minfasta = 200000, output_tag=None, output_compression='none'):
    '''
    Write binned FASTA files

    Returns: DataFrame with information on the bins
    '''
    from collections import defaultdict
    import pandas as pd

    res = defaultdict(list)
    for label, name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    written = []
    for label, contigs in res.items():
        sizes = [len(contig_seqs[contig]) for contig in contigs]
        whole_bin_bp = sum(sizes)

        if whole_bin_bp >= minfasta:
            if output_tag is None:
                ofname = f'bin.{label}.fa'
            else:
                ofname = f'{output_tag}_{label}.fa'
            ofname = os.path.join(output, ofname)
            if output_compression != 'none':
                ofname += '.' + output_compression
            n50, l50 = n50_l50(sizes)
            written.append([ofname, whole_bin_bp, len(contigs), n50, l50])
            with possibly_compressed_write(ofname) as ofile:
                for contig in contigs:
                    ofile.write(f'>{contig}\n{contig_seqs[contig]}\n')
    return pd.DataFrame(written, columns=['filename', 'nbps', 'nr_contigs', 'N50', 'L50'])



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
    if not contig_dict:
        logger.warning(f'No contigs in {fasta_path}')
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
    envn = env.lower().replace('-', '_')
    known_environments = [
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
            'chicken_caecum',
            'global',
            ]
    if envn in known_environments:
        return os.path.join(os.path.split(__file__)[0], f'{envn}_model.h5')
    else:
        sys.stderr.write(
            f"Error: Expected environment '{env}' does not exist (known ones are {', '.join(known_environments)})\n")
        sys.exit(1)

def concatenate_fasta(fasta_files, min_length, output, separator, output_compression='none'):
    """
    Concatenate multiple FASTA files into one

    Returns name of the concatenated file
    """
    ofname = os.path.join(output, 'concatenated.fa')
    if output_compression != 'none':
        ofname += '.' + output_compression
    with possibly_compressed_write(ofname) as concat_out:
        for fname in fasta_files:
            sample_name = os.path.basename(fname).split('.')[0]
            for h, seq in fasta_iter(fname):
                if separator in h:
                    sys.stderr.write(
                        f"Error in file {fname}: Contig ID '{h}' contains the separator ('{separator}'), please choose another separator.\n")
                    sys.exit(1)
                if len(seq) >= min_length:
                    header = f'{sample_name}{separator}{h}'
                    concat_out.write(f'>{header}\n{seq}\n')
    return ofname

def n50_l50(sizes):
    '''Computes N50 & L50 for a set of contig sizes

    Returns (n50, l50)
    '''
    import numpy as np
    s = np.array(sizes)
    s.sort()
    s = s[::-1]
    l50 = np.where(np.cumsum(s)>=s.sum()//2)[0][0]
    n50 = s[l50]
    return n50, l50+1

