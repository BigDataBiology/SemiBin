from __future__ import annotations
import tempfile
from .utils import maybe_uncompress, fasta_iter, describe_returncode, format_stderr
import os
import sys
import json
import subprocess
import shutil
from collections import defaultdict
from typing import IO

_MARKER_CACHE_SCHEMA_VERSION = 1


def _marker_cache_fingerprint(
        fasta_path: str,
        binned_length: int,
        orf_finder: str,
        prodigal_output_faa: str | None) -> dict:
    '''Build a fingerprint of inputs that affect ``markers.hmmout``.

    Uses (size, mtime_ns) of the FASTA rather than a content hash to keep the
    check O(1) on multi-GB files.
    '''
    st = os.stat(fasta_path)
    return {
        'schema_version': _MARKER_CACHE_SCHEMA_VERSION,
        'fasta_abspath': os.path.abspath(fasta_path),
        'fasta_size': st.st_size,
        'fasta_mtime_ns': st.st_mtime_ns,
        'binned_length': int(binned_length),
        'orf_finder': orf_finder,
        'prodigal_output_faa': prodigal_output_faa,
    }

__normalize_marker_trans = {
    'TIGR00388': 'TIGR00389',
    'TIGR00471': 'TIGR00472',
    'TIGR00408': 'TIGR00409',
    'TIGR02386': 'TIGR02387',
}

MARKER_HMM_PATH = os.path.join(os.path.dirname(__file__), 'marker.hmm')

# Number of distinct single-copy marker genes in MARKER_HMM_PATH, after the
# normalization above merges paired TIGRFAM profiles (111 profiles -> 107 genes).
# Used as the denominator when estimating bin completeness.
N_SINGLE_COPY_MARKERS = 107

def get_marker(hmmout: str | IO[str],
               fasta_path: str | None = None,
               min_contig_len: int | None = None,
               multi_mode: bool = False,
               orf_finder: str | None = None,
               contig_to_marker: bool = False,
               ) -> list[str] | dict[str, list[str]]:
    '''Parse HMM output file and return markers

    Hits covering less than 40% of the HMM profile are discarded.

    Parameters
    ----------
    hmmout : str or file-like
        Output of ``hmmsearch --domtblout``
    fasta_path : str, optional
        Contig FASTA file (needed only if ``min_contig_len`` is given)
    min_contig_len : int, optional
        If given, hits on contigs shorter than this are discarded
    multi_mode : bool, optional
        If True, contig names are expected to be of the form
        ``<bin>.<contig>`` and seeds are computed per bin
    orf_finder : str, optional
        ORF finder used to generate the ORFs (determines how ORF names are
        mapped back to contig names)
    contig_to_marker : bool, optional
        If True, return the mapping of contigs to marker genes instead of seeds

    Returns
    -------
    seeds : list of str or dict
        - if ``contig_to_marker``: contig -> list of marker genes (a
          ``defaultdict(list)``)
        - if ``multi_mode``: bin -> list of seed contigs
        - otherwise: list of seed contigs (each containing a copy of the
          same marker gene, so they are presumably from different genomes)

    If no marker hits survive filtering, returns an empty result of the
    appropriate type for the mode (empty mapping for ``contig_to_marker``
    and ``multi_mode``, empty list otherwise)
    '''
    import pandas as pd
    def empty_result():
        if contig_to_marker:
            return defaultdict(list)
        return {} if multi_mode else []

    data = pd.read_table(hmmout, sep=r'\s+',  comment='#', header=None,
                         usecols=(0,3,5,15,16), names=['orf', 'gene', 'qlen', 'qstart', 'qend'])
    if not len(data):
        return empty_result()
    data['gene'] = data['gene'].map(lambda m: __normalize_marker_trans.get(m , m))
    qlen = data[['gene','qlen']].drop_duplicates().set_index('gene')['qlen']

    def contig_name(orf_name: str) -> str:
        if orf_finder in ['prodigal', 'fast-naive']:
            contig,_ = orf_name.rsplit( '_', 1)
        else:
            contig,_,_,_ = orf_name.rsplit( '_', 3)
        return contig

    data = data.query('(qend - qstart) / qlen > 0.4').copy()
    data['contig'] = data['orf'].map(contig_name)
    if min_contig_len is not None:
        contig_len = {h:len(seq) for h,seq in fasta_iter(fasta_path)}
        data = data[data['contig'].map(lambda c: contig_len[c] >= min_contig_len)]
    data = data.drop_duplicates(['gene', 'contig'])
    if not len(data):
        return empty_result()

    if contig_to_marker:
        markers = data['gene'].values
        contigs = data['contig'].values
        contig2markers = defaultdict(list)
        for m, c in zip(markers, contigs, strict=True):
            contig2markers[c].append(m)
        return contig2markers
    else:
        def extract_seeds(gene_counts, hits):
            '''Pick the contigs containing the marker gene with median count

            gene_counts : pd.Series
                gene -> number of contigs with a hit to that gene
            hits : pd.DataFrame
                marker hits (with ``gene`` and ``contig`` columns)
            '''
            gene_counts = gene_counts.sort_values()
            median = gene_counts.iloc[len(gene_counts) //2]

            # the original version broke ties by picking the shortest query, so we
            # replicate that here:
            candidates = gene_counts.index[gene_counts == median]
            seed_gene = qlen.loc[candidates].idxmin()
            seeds = list(hits.query('gene == @seed_gene')['contig'])
            seeds.sort()
            return seeds


        if multi_mode:
            data['bin'] = data['orf'].str.split(pat='.', n=0, expand=True)[0]
            counts = data.groupby(['bin', 'gene'])['orf'].count()
            bin_seeds = {}
            for b, bin_counts in counts.groupby(level=0):
                seeds = extract_seeds(bin_counts.droplevel(0), data.query('bin == @b', local_dict={'b':b}))
                bin_seeds[b] = [c.split('.',1)[1] for c in seeds]
            return bin_seeds
        else:
            counts = data.groupby('gene')['orf'].count()
            return extract_seeds(counts, data)

def estimate_seeds(fasta_path: str,
                 binned_length: int,
                 num_process: int,
                 multi_mode: bool = False,
                 output: str | None = None,
                 orf_finder: str = 'prodigal',
                 prodigal_output_faa: str | None = None,
                 ) -> list[str] | dict[str, list[str]]:
    '''Estimate seed contigs (and thus, the number of bins) from a FASTA file

    Runs the ORF finder and ``hmmsearch`` against single-copy marker genes;
    the contigs containing a copy of the marker gene with the median count
    are used as seeds (see ``get_marker``).

    Parameters
    ----------
    fasta_path : str
        Input contigs (FASTA, possibly compressed)
    binned_length : int
        Minimal contig length
    num_process : int
        Number of CPUs to use
    multi_mode : bool, optional
        If True, treat input as resulting from concatenating multiple bins
        (contigs named ``<bin>.<contig>``) and return seeds per bin
    output : str, optional
        If given, ``markers.hmmout`` is saved to this directory (and reused on
        later calls if the inputs have not changed)
    orf_finder : str, optional
        One of 'prodigal', 'fast-naive', or 'fraggenescan'
    prodigal_output_faa : str, optional
        Pre-computed ORFs to use instead of running the ORF finder

    Returns
    -------
    seeds : list of str or dict
        See ``get_marker``
    '''
    import logging
    from .orffinding import run_orffinder
    from .atomicwrite import atomic_write
    logger = logging.getLogger("SemiBin2")
    with tempfile.TemporaryDirectory() as tdir:
        fingerprint = _marker_cache_fingerprint(
                fasta_path, binned_length, orf_finder, prodigal_output_faa)
        input_fasta = fasta_path
        fasta_path = maybe_uncompress(fasta_path, tdir)
        if output is not None:
            os.makedirs(output, exist_ok=True)
            target_dir = output
            hmm_output = os.path.join(output, 'markers.hmmout')
            sidecar = hmm_output + '.json'
            if os.path.exists(hmm_output):
                cached = None
                if os.path.exists(sidecar):
                    try:
                        with open(sidecar) as f:
                            cached = json.load(f)
                    except (OSError, ValueError) as e:
                        logger.warning(
                            f'Could not read marker cache metadata at {sidecar} ({e}); recomputing.')
                if cached == fingerprint:
                    return get_marker(hmm_output,
                                        fasta_path, binned_length, multi_mode, orf_finder=orf_finder)
                if cached is None:
                    logger.warning(
                        f'Found {hmm_output} without matching cache metadata; recomputing markers.')
                else:
                    logger.warning(
                        f'Inputs to marker calling changed since {hmm_output} was written; recomputing markers.')
                os.remove(hmm_output)
                if os.path.exists(sidecar):
                    os.remove(sidecar)
        else:
            target_dir = tdir

        # Check before running the ORF finder (which can take a while)
        hmmsearch_path = shutil.which('hmmsearch')
        if hmmsearch_path is None:
            logger.error('hmmsearch binary not found in PATH (it is required to identify marker genes). '
                         'Please install HMMER.')
            sys.exit(1)

        contig_output = run_orffinder(fasta_path, num_process, tdir, orf_finder, prodigal_output_faa=prodigal_output_faa)

        hmm_output = os.path.join(target_dir, 'markers.hmmout')
        try:
            with open(os.path.join(tdir, 'markers.hmmout.out'), 'w') as hmm_out_log, \
                    open(os.path.join(tdir, 'markers.hmmout.err'), 'w') as hmm_err_log:
                subprocess.check_call(
                    [hmmsearch_path,
                     '--domtblout',
                     hmm_output,
                     '--cut_tc',
                     '--cpu', str(num_process),
                     MARKER_HMM_PATH,
                     contig_output,
                     ],
                    stdout=hmm_out_log,
                    stderr=hmm_err_log,
                )
        except (OSError, subprocess.CalledProcessError) as e:
            if os.path.exists(hmm_output):
                os.remove(hmm_output)
            if isinstance(e, subprocess.CalledProcessError):
                with open(os.path.join(tdir, 'markers.hmmout.err')) as hmm_err_log:
                    err = format_stderr(hmm_err_log.read())
                logger.error(f'hmmsearch {describe_returncode(e.returncode)} '
                             f'(when searching for marker genes in ORFs predicted from `{input_fasta}`).{err}')
            else:
                logger.error(f'Running hmmsearch failed: {e}')
            sys.exit(1)

        if output is not None:
            with atomic_write(hmm_output + '.json', mode='wt', overwrite=True) as f:
                json.dump(fingerprint, f, sort_keys=True)

        return get_marker(hmm_output, fasta_path, binned_length, multi_mode, orf_finder=orf_finder)

