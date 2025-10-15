import os
import subprocess
import multiprocessing
import sys
import random
import contextlib
import pathlib

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
    args.training_type = 'auto'

    if args.training_type == 'semi' and args.self_supervised:
        logger.error('Both --training-type=semi and --self-supervised arguments used')
        sys.exit(1)

    if args.training_type == 'self' and args.semi_supervised:
        logger.warning('Both --training-type=self and --semi-supervised arguments used')
        sys.exit(1)

    if not args.self_supervised and not args.semi_supervised:
        logger.debug(
            f"SemiBin will run in self supervised mode")
        args.training_type = 'self'

    elif args.self_supervised and args.semi_supervised:
        logger.warning(
            f'You chose both semi-supervised and self-supervised learning! SemiBin will use self-supervised learning')
        args.training_type = 'self'

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
                logger.error(
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
        os.environ['OPENBLAS_NUM_THREADS'] = str(args.num_process)
        os.environ['MKL_NUM_THREADS'] = str(args.num_process)
        os.environ['VECLIB_MAXIMUM_THREADS'] = str(args.num_process)

    if args.cmd in ['train', 'train_semi']:
        args.cmd = 'train_semi'
        args.training_type = 'semi'

    if args.cmd == 'train_self':
        args.training_type = 'self'

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train_semi', 'bin']:
        if args.orf_finder not in ['prodigal', 'fraggenescan', 'fast-naive']:
            logger.error(
                f"Error: SemiBin only supports 'prodigal'/'fraggenescan'/'fast-naive' as the ORF finder (--orf-finder option).\n")
            exit_with_error = True
        if args.orf_finder == 'fraggenescan':
            from time import sleep
            logger.warning('Using FragGeneScan (--orf-finder=fraggenescan) as an ORF finder is deprecated and will be removed in the next version of SemiBin. '
                        'Please use fast-naive (the default) or prodigal instead.')
            if sys.stdout.isatty():
                for i in range(5, 0, -1):
                    print(f'SemiBin will continue in {i} seconds...', end='\r')
                    sleep(1)
        expect_file(args.prodigal_output_faa)
        if args.prodigal_output_faa is not None:
            logger.warning('Using --prodigal-output-faa is deprecated and will be removed in a future version of SemiBin')
            args.orf_finder = 'prodigal'

    if hasattr(args, 'engine'):
        args.engine = args.engine.lower()
        if args.engine not in ['auto', 'gpu', 'cpu']:
            logger.error(
                f"Error: Argument '--engine' needs to be one of auto[default]/gpu/cpu.\n")
            exit_with_error = True

    if args.cmd == 'generate_cannot_links':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.taxonomy_results_fname)

    if args.cmd in ['generate_sequence_features_single', 'generate_sequence_features_multi', 'single_easy_bin', 'multi_easy_bin']:
        if args.bams and args.abundances:
            logger.error(
                f"Error: can not use BAM files and abundance files at the same time.\n")
            exit_with_error = True

    if args.cmd == 'generate_sequence_features_single':
        expect_file(args.contig_fasta)
        if args.bams:
            expect_file_list(args.bams)
        if args.abundances:
            expect_file_list(args.abundances)

    if args.cmd == 'generate_sequence_features_multi':
        expect_file(args.contig_fasta)
        if args.bams:
            expect_file_list(args.bams)
        if args.abundances:
            expect_file_list(args.abundances)

    if args.cmd in ['train_semi', 'train_self']:
        if not args.train_from_many:
            if len(args.data) > 1:
                logger.error(
                    f"Error: Expected one data.csv file with single mode.\n")
                exit_with_error = True

            if len(args.data_split) > 1:
                logger.error(
                    f"Error: Expected one data_split.csv file with single mode.\n")
                exit_with_error = True
            if args.cmd == 'train_semi':
                if len(args.contig_fasta) > 1:
                    logger.error(
                        f"Error: Expected one fasta file with single mode.\n")
                    exit_with_error = True

                if len(args.cannot_link) > 1:
                    logger_error(
                        f"Error: Expected one cannot.txt file with single mode.\n")
                    exit_with_error = True
                expect_file(args.cannot_link[0])
                expect_file(args.contig_fasta[0])

            expect_file(args.data[0])
            expect_file(args.data_split[0])

        else:
            if args.cmd == 'train_semi':
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
            logger.error(
                f"Error: Did not understand sequencing_type argument '{args.sequencing_type}' (should be short_reads or long_reads).\n")
            exit_with_error = True
        logger.info(f'Binning for {args.sequencing_type}')

    if args.cmd == 'bin':
        if args.environment is None and args.model_path is None:
            logger.error(
                f"Error: Please choose input a model path or use our built-in model for [human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global].\n")
            exit_with_error = True
        if args.environment is not None and args.model_path is not None:
            logger.error(
                f"Error: You cannot use both an explicit model path and an environment.\n")
            exit_with_error = True
        if args.model_path is not None:
            expect_file(args.model_path)
        expect_file(args.contig_fasta)
        expect_file(args.data)

    if getattr(args, 'environment', None) is not None:
        try:
            get_model_path(args.environment)
        except KeyError as e:
            logger.error(e.args[0])
            exit_with_error = True
    if args.cmd == 'single_easy_bin':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        if args.bams:
            expect_file_list(args.bams)
        if args.abundances:
            expect_file_list(args.abundances)

        else:
            check_training_type(logger, args)

    if args.cmd == 'multi_easy_bin':
        check_training_type(logger, args)
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        if args.bams:
            expect_file_list(args.bams)
        if args.abundances:
            expect_file_list(args.abundances)

        if args.training_type not in ['semi', 'self']:
            logger.error(
                f"Error: You need to specify the training algorithm in semi/self.\n")
            exit_with_error = True

    if args.cmd == 'concatenate_fasta':
        expect_file_list(args.contig_fasta)
        contig_name = []
        for contig in args.contig_fasta:
            contig_name.append(os.path.basename(contig).split('.')[0])
        if len(set(contig_name)) != len(contig_name):
            logger.error(
                f"Error: Make sure that every contig file have different names.\n")
            exit_with_error = True


    if getattr(args, 'train_from_many', False):
        args.mode = 'several'
    elif args.cmd in ['train_semi', 'train_self'] and not hasattr(args, 'mode'):
        args.mode = 'single'

    if getattr(args, 'write_pre_reclustering_bins', False) and \
            not getattr(args, 'recluster', True):
        logger.error(
            f"Error: Cannot use --write-pre-reclustering-bins with --no-recluster.\n")
        exit_with_error = True

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



def maybe_uncompress(fafile, tdir):
    if fafile.endswith('.gz') or \
            fafile.endswith('.bz2') or \
            fafile.endswith('.xz'):
        oname = f'{tdir}/expanded.fa'
        with open(oname, 'wt') as out:
            for header, seq in fasta_iter(fafile):
                out.write(f'>{header}\n{seq}\n')
        return oname
    return fafile



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

def load_fasta(fasta_path: str, ratio: float):
    """
    Returns

    computed_min_length: minimum length of contigs (1000 or 2500, depending on the ratio)
    must_link_threshold: threshold to break up contigs
    contigs: dictionary ID -> contig sequence
    """
    total_bps = 0
    contig_bp_2500 = 0
    contig_length_list = []
    contig_dict = {}

    for h, seq in fasta_iter(fasta_path):
        if 1000 <= len(seq) <= 2500:
            contig_bp_2500 += len(seq)
        contig_length_list.append(len(seq))
        total_bps += len(seq)
        contig_dict[h] = seq

    computed_min_length = (
                1000
                if contig_bp_2500 / total_bps < ratio
                else 2500)
    must_link_threshold = get_must_link_threshold(contig_length_list)
    if not contig_dict:
        import logging
        logger = logging.getLogger('SemiBin2')
        logger.warning(f'No contigs in {fasta_path}')
    return computed_min_length, must_link_threshold, contig_dict


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

def get_model_path(env : str) -> str:
    '''Finds model path

    If not found, triggers `sys.exit(1)`
    '''
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
        # From Python 3.9, we can use importlib.resources
        return os.path.join(os.path.split(__file__)[0], 'models', f'{envn}_model.pt')
    else:
        raise KeyError(
                f"Error: Unknown environment '{env}' does not exist (known ones are {', '.join(known_environments)})")

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


def maybe_crams2bams(bams, contig_fasta : str, num_process : int, odir : str): # bams : list[str] is not available on Python 3.7
    '''
    maybe_crams2bams converts CRAM to BAM

    Parameters
    ----------
    bams : list of str
        List of BAM/CRAM files
    contig_fasta : str
        Contig FASTA file
    num_process : int
        Number of processes to use
    odir : str
        Output directory for extracted BAM files

    Returns
    -------
    obams : list of str
        List of BAM files (extracted if CRAM or original)
    '''
    rs = []
    for bam in bams:
        if bam.endswith('.cram'):
            obam = pathlib.PurePath(odir) / pathlib.PurePath(bam).with_suffix('.bam').name
            with open(obam, 'wb') as cram_out:
                subprocess.check_call(
                    ['samtools', 'view',
                     '-bS',
                     '-@', str(num_process),
                     '-T', contig_fasta,
                     bam],
                    stdout=cram_out)
            rs.append(str(obam))
        else:
            rs.append(bam)
    return rs


def maybe_compute_min_length(min_length, fafile, ratio):
    if min_length is not None:
        return min_length
    c_min_len, _, _ = load_fasta(fafile, ratio)
    return c_min_len


def norm_abundance(data, features):
    import pandas as pd
    import numpy as np
    assert isinstance(data, pd.DataFrame), "norm_abundance did not receive a pandas dataframe"
    n = data.shape[1] - len(features["kmer"]) - len(features["motif"]) - len(features["motif_present"])
    assert n == len(features["depth"]), "Depth should equal all_features - motifs - kmer!"
    flag = False

    if n >= 20:
        flag = True
    else:
        flag = (n >= 5 and np.mean(np.sum(data[features["depth"]], axis=1)) > 2)

    return flag



    
import re

# IUPAC codes dictionary for sequence pattern matching
iupac_chars = [
    "A", "T", "C", "G",
    "R", "Y", 
    "S", "W", 
    "K", "M",
    "B",
    "D",
    "H",
    "V",
    "N"
]

def check_motif(column_name):
    """
    Check if a column is a motif.

    Parameters:
    column_name (str): The column name to check.

    Returns:
    bool: True if the column is a motif, False otherwise.
    """
    # Regex pattern to capture the motif
    pattern = r"(\b\w+_(m|a|21839)-\d+)"
    
    motif_mod = re.search(pattern, str(column_name))
   
    if not motif_mod:
        return False

    motif_mod = motif_mod.group(1)
    
    string_motif, mod_pos = motif_mod.rsplit('_', 1)  # rsplit on the last underscore
    _, motif = string_motif.rsplit('_', 1)
    mod, pos = mod_pos.split('-')

    if all(char in iupac_chars for char in motif):
        # Check if mod is valid and pos is within range
        if mod in ["m", "a", "21839"] and 0 <= int(pos) < 20:
            return True
        else:
            return False
    else:
        return False
    

def get_features(df):
    """
    Takes a DataFrame and extracts indices of features to populate the provided features dictionary.
    Specific features are extracted based on the column names:
    - Indices of columns ending in 'bam_mean' or 'bam_var' are considered depth features.
    
    Parameters:
        df (pd.DataFrame): The DataFrame from which to extract features.
        features_dict (dict): The dictionary to populate with feature indices.
    """
    features_dict = {
        'kmer': [str(i) for i in range(136)],
        'depth': [],
        'motif': [],
        'motif_present': []
    }
    
    columns = df.columns
    # Populate 'depth' with indices of columns ending with 'bam_mean' or 'bam_var'
    features_dict['depth'] = [column for i, column in enumerate(columns) if column.endswith('mean') or column.endswith('var') or column.endswith('bam_cov') or column.endswith(".txt")]
    
    # Populate 'motif' with indices of columns
    features_dict['motif'] = [column for i, column in enumerate(columns) if column.startswith("methylation_value") and check_motif(column)]
    
    features_dict['motif_present'] = [column for i, column in enumerate(columns) if column.startswith("motif_present") and check_motif(column)]
    assert len(features_dict['depth'] + features_dict['motif'] + features_dict['motif_present']) == len(set(features_dict['depth'] + features_dict['motif'] + features_dict['motif_present']) )
    return features_dict

def normalize_kmer_motif_features(train_data, train_data_split):
    """
    MinMax scaling was chosen to normalize kmers which is in the range of 0-0.07 however motifs are in the range of 0-1.
    """
    import numpy as np
    from sklearn.preprocessing import MinMaxScaler
   
    scaler = MinMaxScaler()
    train_data = scaler.fit_transform(train_data)
    train_data_split = scaler.transform(train_data_split)

    assert train_data.shape[1] == train_data_split.shape[1]
    # assert (train_data.shape[1] - train_data_split.shape[1]) % 2 == 0
    return train_data, train_data_split
