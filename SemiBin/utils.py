from __future__ import annotations
import os
import subprocess
import multiprocessing
import logging
import sys
import random
import contextlib
import pathlib
import argparse
from collections.abc import Iterator, Sequence
from typing import TYPE_CHECKING, TextIO

from .fasta import fasta_iter

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd


@contextlib.contextmanager
def possibly_compressed_write(filename: str) -> Iterator[TextIO]:
    '''Atomically write a text file, compressing it based on its extension

    Compression is inferred from the extension (``.gz``, ``.bz2``, or ``.xz``;
    anything else is written uncompressed). The file is only moved into place
    if the block exits without an exception.

    Usage::

        with possibly_compressed_write('output.fa.gz') as out:
            out.write('>contig\nACGT\n')
    '''
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

def describe_returncode(rc: int) -> str:
    '''Human-readable description of a (failed) subprocess return code'''
    if rc < 0:
        import signal
        try:
            sig = signal.Signals(-rc).name
        except ValueError:
            sig = f'signal {-rc}'
        return f'was killed by {sig}'
    return f'exited with return code {rc}'


def format_stderr(data: str | bytes, max_lines: int = 10) -> str:
    '''Format the tail of a captured stderr for inclusion in an error message

    Parameters
    ----------
    data : str or bytes
        Captured stderr
    max_lines : int, optional
        Only the last ``max_lines`` lines are included

    Returns
    -------
    formatted : str
        Text to append to an error message (an empty string if there was no
        output)
    '''
    if isinstance(data, bytes):
        data = data.decode('utf-8', errors='replace')
    lines = data.strip().splitlines()
    if not lines:
        return ''
    if len(lines) > max_lines:
        lines = ['[...]'] + lines[-max_lines:]
    return '\nError output (stderr):\n' + '\n'.join('    ' + ell for ell in lines)


def set_training_type(logger: logging.Logger, args: argparse.Namespace) -> None:
    '''Set ``args.training_type`` to 'self' or 'semi' based on the command line flags

    Self-supervised learning is the default (and wins if both are requested).
    '''
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


def validate_normalize_args(logger: logging.Logger, args: argparse.Namespace) -> None:
    '''Validate and normalize command line arguments

    ``args`` is modified in place (e.g., aliases are resolved and defaults
    that depend on other arguments are filled in). If any errors are found,
    they are all logged and then the program exits with ``sys.exit(1)``.
    '''
    exit_with_error = False

    def expect_file(f: str | None) -> None:
        nonlocal exit_with_error
        if f is not None:
            if not os.path.exists(f):
                logger.error(
                    f"Error: Expected file '{f}' does not exist\n")
                exit_with_error = True

    def expect_file_list(fs: list[str] | None) -> None:
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

    if args.cmd in ['single_easy_bin', 'multi_easy_bin', 'train_semi', 'bin', 'bin_long']:
        if args.orf_finder not in ['prodigal', 'fraggenescan', 'fast-naive']:
            logger.error(
                "SemiBin only supports 'prodigal'/'fraggenescan'/'fast-naive' as the ORF finder (--orf-finder option).")
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
                "Argument '--engine' needs to be one of auto[default]/gpu/cpu.")
            exit_with_error = True

    if args.cmd == 'generate_cannot_links':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.taxonomy_results_fname)

    if args.cmd in ['generate_sequence_features_single', 'generate_sequence_features_multi', 'single_easy_bin', 'multi_easy_bin']:
        if args.bams and args.abundances:
            logger.error(
                "Cannot use BAM files and abundance files at the same time.")
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
                    "Expected one data.csv file with single mode.")
                exit_with_error = True

            if len(args.data_split) > 1:
                logger.error(
                    "Expected one data_split.csv file with single mode.")
                exit_with_error = True
            if args.cmd == 'train_semi':
                if len(args.contig_fasta) > 1:
                    logger.error(
                        "Expected one fasta file with single mode.")
                    exit_with_error = True

                if len(args.cannot_link) > 1:
                    logger.error(
                        "Expected one cannot.txt file with single mode.")
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
        elif args.sequencing_type.lower() in ['long', 'long-read', 'long_reads', 'long-reads', 'long_read', 'hybrid']:
            args.sequencing_type = 'long_read'
        else:
            logger.error(
                f"Did not understand sequencing_type argument '{args.sequencing_type}' (should be short_reads or long_reads).")
            exit_with_error = True
        logger.info(f'Binning for {args.sequencing_type}')

    if args.cmd in ['bin', 'bin_long']:
        if args.environment is None and args.model_path is None:
            logger.error(
                "Please provide a model path or use one of our built-in models: [human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global].")
            exit_with_error = True
        if args.environment is not None and args.model_path is not None:
            logger.error(
                "You cannot use both an explicit model path and an environment.")
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
        set_training_type(logger, args)
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        if args.bams:
            expect_file_list(args.bams)
        if args.abundances:
            expect_file_list(args.abundances)

    if args.cmd == 'multi_easy_bin':
        set_training_type(logger, args)
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        if args.bams:
            expect_file_list(args.bams)
        if args.abundances:
            expect_file_list(args.abundances)

        if args.training_type not in ['semi', 'self']:
            logger.error(
                "You need to specify the training algorithm as semi/self.")
            exit_with_error = True

    if args.cmd == 'concatenate_fasta':
        expect_file_list(args.contig_fasta)
        sample_names = []
        for fasta_file in args.contig_fasta:
            sample_names.append(os.path.basename(fasta_file).split('.')[0])
        if len(set(sample_names)) != len(sample_names):
            logger.error(
                "Make sure that every contig file has a different name.")
            exit_with_error = True


    if getattr(args, 'train_from_many', False):
        args.mode = 'several'
    elif args.cmd in ['train_semi', 'train_self'] and not hasattr(args, 'mode'):
        args.mode = 'single'

    if getattr(args, 'write_pre_reclustering_bins', False) and \
            not getattr(args, 'recluster', True):
        logger.error(
            "Cannot use --write-pre-reclustering-bins with --no-recluster.")
        exit_with_error = True

    if exit_with_error:
        sys.exit(1)


def get_must_link_threshold(contig_lens: Sequence[int]) -> int:
    """
    Calculate the minimum length for contigs to be split into must-link pairs

    The threshold is the length such that contigs shorter than it account for
    approximately 2% of the total assembly size, but never less than 4,000bp.
    """
    import numpy as np
    contig_lens = np.array(contig_lens)
    contig_lens.sort()
    frac = np.cumsum(contig_lens)/np.sum(contig_lens)
    ix = np.argmax(frac > 0.02) # argmax finds first True element
    threshold = contig_lens[ix]
    return np.clip(threshold, 4000, None)

def parse_mmseqs(mmseqs_result: pd.DataFrame
                 ) -> tuple[list[tuple[str, str]], list[tuple[str, str]], list[tuple[str, str]]]:
    '''Generate cannot-link pairs from taxonomic annotations

    Two contigs confidently annotated to different species (or different
    genera) cannot be in the same bin.

    Parameters
    ----------
    mmseqs_result : pd.DataFrame
        Output of ``mmseqs taxonomy`` with (at least) the columns
        ``contig_name``, ``rank_name``, ``scientific_name``, ``score``, and
        ``lineage``

    Returns
    -------
    cannot_link_species : list of (str, str)
        Pairs of contigs annotated to different species (score > 0.95)
    cannot_link_genus : list of (str, str)
        Pairs of contigs annotated to different genera (score > 0.80)
    cannot_link_mix : list of (str, str)
        Pairs (species-annotated contig, genus-annotated contig) where the
        genus of the first does not match the second
    '''
    species_hits = mmseqs_result.query('rank_name == "species" and score > 0.95')
    genus_hits = mmseqs_result.query('rank_name == "genus" and score > 0.80')
    species_contigs = species_hits['contig_name'].tolist()
    species_names = species_hits['scientific_name'].tolist()
    species_lineages = species_hits['lineage'].tolist()
    genus_contigs = genus_hits['contig_name'].tolist()
    genus_names = genus_hits['scientific_name'].tolist()

    cannot_link_species = []
    for i in range(len(species_contigs)):
        for j in range(i + 1, len(species_contigs)):
            if species_names[i] != species_names[j]:
                cannot_link_species.append(
                    (species_contigs[i], species_contigs[j]))

    cannot_link_genus = []
    for i in range(len(genus_contigs)):
        for j in range(i + 1, len(genus_contigs)):
            if genus_names[i] != genus_names[j]:
                cannot_link_genus.append(
                    (genus_contigs[i], genus_contigs[j]))

    cannot_link_mix = []
    for species_contig, lineage in zip(species_contigs, species_lineages, strict=True):
        genus_of_species = lineage.split(';')[-2]
        for genus_contig, genus_name in zip(genus_contigs, genus_names, strict=True):
            if genus_of_species != genus_name:
                cannot_link_mix.append(
                    (species_contig, genus_contig))

    return cannot_link_species, cannot_link_genus, cannot_link_mix


def generate_cannot_link(
        mmseqs_path: str,
        namelist: list[str],
        n_must_link: int,
        output: str,
        sample: str) -> None:
    '''Generate cannot-link constraints from taxonomic annotations

    Writes ``<output>/<sample>.txt`` with one comma-separated pair of contig
    names per line. The number of constraints is subsampled to be
    proportional to the number of must-link pairs (capped at 4 million).

    Parameters
    ----------
    mmseqs_path : str
        Taxonomy results (``mmseqs createtsv`` output)
    namelist : list of str
        Contigs to consider (others are ignored)
    n_must_link : int
        Number of contigs that will be split into must-link pairs
    output : str
        Output directory
    sample : str
        Sample name (output basename)
    '''
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
    num_whole_data = np.clip(1000 * n_must_link, None, 4_000_000)
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



def maybe_uncompress(fafile: str, tdir: str) -> str:
    '''Uncompress a FASTA file into ``tdir`` if it is compressed

    Returns the path to an uncompressed FASTA file (``fafile`` itself if it
    was not compressed).
    '''
    if fafile.endswith('.gz') or \
            fafile.endswith('.bz2') or \
            fafile.endswith('.xz'):
        oname = f'{tdir}/expanded.fa'
        with open(oname, 'wt') as out:
            for header, seq in fasta_iter(fafile):
                out.write(f'>{header}\n{seq}\n')
        return oname
    return fafile



def write_bins(namelist: list[str],
               contig_labels: Sequence[int],
               output: str,
               contig_seqs: dict[str, str],
               minfasta: int = 200000,
               output_tag: str | None = None,
               output_compression: str = 'none') -> pd.DataFrame:
    '''
    Write binned FASTA files

    Parameters
    ----------
    namelist : list of str
        Contig names
    contig_labels : sequence of int
        Bin label for each contig in ``namelist`` (-1 means unbinned)
    output : str
        Output directory
    contig_seqs : dict[str, str]
        contig name -> sequence
    minfasta : int, optional
        Bins whose total size (in basepairs) is smaller than this are not
        written
    output_tag : str, optional
        If given, bins are named ``<output_tag>_<label>.fa``, otherwise
        ``bin.<label>.fa``
    output_compression : str, optional
        'none' or one of 'gz', 'bz2', 'xz'

    Returns
    -------
    bins_info : pd.DataFrame
        Information on the written bins (columns: ``filename``, ``nbps``,
        ``nr_contigs``, ``N50``, ``L50``)
    '''
    from collections import defaultdict
    import pandas as pd

    bin_contigs = defaultdict(list)
    for label, name in zip(contig_labels, namelist, strict=True):
        if label != -1:
            bin_contigs[label].append(name)

    os.makedirs(output, exist_ok=True)

    written = []
    for label, contigs in bin_contigs.items():
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


def log_binning_stats(
        logger: logging.Logger,
        bins_info: pd.DataFrame,
        contig_dict: dict[str, str],
        binned_length: int | None = None) -> None:
    '''Log summary statistics about the produced bins

    Parameters
    ----------
    logger
        logger to write the summary to
    bins_info : pd.DataFrame
        the DataFrame returned by ``write_bins`` (one row per bin, with
        ``nbps`` and ``nr_contigs`` columns)
    contig_dict : dict
        contig-id -> sequence for the whole input assembly
    binned_length : int, optional
        minimum contig length considered for binning; when given, contigs
        shorter than this are reported as a separate "too short" category
    '''
    n_bins = len(bins_info)
    if n_bins == 0:
        return

    total_contigs = len(contig_dict)
    total_bp = sum(len(s) for s in contig_dict.values())

    binned_contigs = int(bins_info['nr_contigs'].sum())
    binned_bp = int(bins_info['nbps'].sum())

    def _pct(part, whole):
        return 100.0 * part / whole if whole else 0.0

    if binned_length is not None:
        short_contigs = sum(len(s) < binned_length for s in contig_dict.values())
        short_bp = sum(len(s) for s in contig_dict.values() if len(s) < binned_length)
    else:
        short_contigs = short_bp = 0

    logger.info('Binning summary:')
    logger.info(f'  Bins produced: {n_bins}')
    logger.info(
        f'  Contigs binned: {binned_contigs:,} / {total_contigs:,} '
        f'({_pct(binned_contigs, total_contigs):.1f}%)')
    logger.info(
        f'  Basepairs binned: {binned_bp:,} / {total_bp:,} '
        f'({_pct(binned_bp, total_bp):.1f}%)')
    if binned_length is not None:
        eligible_unbinned = total_contigs - binned_contigs - short_contigs
        logger.info(
            f'  Unbinned: {eligible_unbinned:,} contigs '
            f'({total_bp - binned_bp - short_bp:,} bp) passed the length '
            f'filter but were not binned; '
            f'{short_contigs:,} contigs ({short_bp:,} bp) too short '
            f'(< {binned_length:,} bp)')
    else:
        logger.info(
            f'  Unbinned: {total_contigs - binned_contigs:,} contigs '
            f'({total_bp - binned_bp:,} bp)')
    sizes = bins_info['nbps']
    logger.info(
        f'  Bin size (bp): mean {binned_bp // n_bins:,}, '
        f'median {int(sizes.median()):,}, '
        f'min {int(sizes.min()):,}, max {int(sizes.max()):,}')
    logger.info(
        f'  Contigs per bin: mean {binned_contigs / n_bins:.1f}, '
        f'median {int(bins_info["nr_contigs"].median())}')
    logger.info(
        '  Note: these bins are not quality-controlled and should not be '
        'treated as MAGs without further QC (e.g., checking completeness and '
        'contamination with CheckM2, and chimerism with GUNC)')


def set_random_seed(seed: int) -> None:
    '''Seed all random number generators used (Python, NumPy, and PyTorch)'''
    import torch
    import numpy as np
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    np.random.seed(seed)

def load_fasta(fasta_path: str, ratio: float) -> tuple[int, int, dict[str, str]]:
    """Load contigs and compute length thresholds

    Parameters
    ----------
    fasta_path : str
        Contig FASTA file
    ratio : float
        If contigs of 1000-2500bp make up less than this fraction of the
        assembly (in basepairs), the minimum length is 1000bp, otherwise 2500bp

    Returns
    -------
    computed_min_length : int
        minimum length of contigs (1000 or 2500, depending on the ratio)
    must_link_threshold : int
        threshold to break up contigs (see ``get_must_link_threshold``)
    contigs : dict[str, str]
        dictionary ID -> contig sequence

    If the file is empty, returns default thresholds and an empty dictionary.
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

    if total_bps == 0:
        logger = logging.getLogger('SemiBin2')
        logger.warning(f'No contigs in {fasta_path}')
        # Caller is expected to detect the empty contig_dict and report a
        # user-facing error; return safe defaults so we do not crash here.
        return 1000, 4000, contig_dict

    computed_min_length = (
                1000
                if contig_bp_2500 / total_bps < ratio
                else 2500)
    must_link_threshold = get_must_link_threshold(contig_length_list)
    return computed_min_length, must_link_threshold, contig_dict


def split_data(data: pd.DataFrame, sample: str, separator: str, is_combined: bool = True) -> pd.DataFrame:
    """
    Extract the rows of a single sample (in multi-sample binning)

    Parameters
    ----------
    data : pd.DataFrame
        Features with a ``contig_name`` column of the form
        ``<sample><separator><contig>``
    sample : str
        Sample to extract
    separator : str
        Separator between sample and contig names
    is_combined : bool, optional
        If True, each column is rescaled (divided by its mean, rounded up to
        a multiple of 100)

    Returns
    -------
    sample_data : pd.DataFrame
        Rows of ``sample``, indexed by contig name (without the sample prefix)
    """
    import numpy as np
    part_data = data[data['contig_name'].str.startswith(
        sample + separator)]
    part_data = part_data.set_index('contig_name')
    part_data.index.name = None
    part_data.index = [ix.split(separator)[1] for ix in part_data.index]
    if is_combined:
        abun_scale = (part_data.mean() / 100).apply(np.ceil) * 100
        part_data = part_data.div(abun_scale)
    return part_data

def get_model_path(env : str) -> str:
    '''Finds the path of the pretrained model for an environment

    Environment names are case-insensitive and dashes are treated as
    underscores (e.g., 'Human-Gut' is the same as 'human_gut').

    Raises ``KeyError`` (with a user-friendly message) if the environment is
    unknown.
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

def concatenate_fasta(fasta_files: list[str],
                      min_length: int,
                      output: str,
                      separator: str,
                      output_compression: str = 'none',
                      keep_full_header: bool = False,
                      ) -> str:
    """
    Concatenate multiple FASTA files into one (for multi-sample binning)

    Contigs are renamed to ``<sample><separator><contig>``, where the sample
    name is the basename of the input file up to the first dot.

    Exits the program if a contig ID already contains ``separator``.

    Parameters
    ----------
    fasta_files : list of str
        Input FASTA files (one per sample)
    min_length : int
        Contigs shorter than this are dropped
    output : str
        Output directory
    separator : str
        Separator between sample and contig names
    output_compression : str, optional
        'none' or one of 'gz', 'bz2', 'xz'
    keep_full_header : bool, optional
        If True, any text in the header after the contig ID is preserved

    Returns
    -------
    ofname : str
        Name of the concatenated file
    """
    logger = logging.getLogger('SemiBin2')
    ofname = os.path.join(output, 'concatenated.fa')
    if output_compression != 'none':
        ofname += '.' + output_compression
    with possibly_compressed_write(ofname) as concat_out:
        for fname in fasta_files:
            sample_name = os.path.basename(fname).split('.')[0]
            for full_header, seq in fasta_iter(fname, full_header=True):
                contig_id, _, description = full_header.partition(' ')
                if separator in contig_id:
                    logger.error(
                        f"In file {fname}: contig ID '{contig_id}' contains the separator ('{separator}'), please choose another separator.")
                    sys.exit(1)
                if len(seq) >= min_length:
                    header = f'{sample_name}{separator}{contig_id}'
                    if keep_full_header and description:
                        header += f' {description}'
                    concat_out.write(f'>{header}\n{seq}\n')
    return ofname

def n50_l50(sizes: Sequence[int]) -> tuple[int, int]:
    '''Computes N50 & L50 for a set of contig sizes

    Returns
    -------
    n50 : int
        Length of the contig such that contigs at least this long cover
        (at least) half the total size
    l50 : int
        Number of contigs at least ``n50`` long (i.e., the minimal number of
        contigs that cover half the total size)
    '''
    import numpy as np
    sorted_sizes = np.array(sizes)
    sorted_sizes.sort()
    sorted_sizes = sorted_sizes[::-1]
    ix50 = np.where(np.cumsum(sorted_sizes) >= sorted_sizes.sum()//2)[0][0]
    n50 = sorted_sizes[ix50]
    return n50, ix50 + 1


def maybe_crams2bams(bams: list[str], contig_fasta: str, num_process: int, odir: str) -> list[str]:
    '''
    maybe_crams2bams converts CRAM to BAM (BAM files are passed through)

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
    obams = []
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
            obams.append(str(obam))
        else:
            obams.append(bam)
    return obams


def maybe_compute_min_length(min_length: int | None, fafile: str, ratio: float) -> int:
    '''Return ``min_length`` if given, otherwise compute it from the contigs

    See ``load_fasta`` for how the minimum length is computed.
    '''
    if min_length is not None:
        return min_length
    c_min_len, _, _ = load_fasta(fafile, ratio)
    return c_min_len


def norm_abundance(data: np.ndarray) -> bool:
    '''Whether abundance features should be normalized

    Parameters
    ----------
    data : np.ndarray
        Feature matrix: 136 k-mer features followed by one abundance column
        per sample

    Returns
    -------
    normalize : bool
        True if there are at least 20 samples, or at least 5 samples with a
        mean total abundance per contig above 2
    '''
    import numpy as np
    n_samples = data.shape[1] - 136

    if n_samples >= 20:
        return True
    return (n_samples >= 5 and np.mean(np.sum(data[:, 136:], axis=1)) > 2)
