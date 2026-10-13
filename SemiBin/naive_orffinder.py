from collections import namedtuple
import multiprocessing as mp
from os import path
from . import fasta

ORFPos = namedtuple('ORFPos', ['start', 'end', 'rc'])

START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']
MIN_LEN = 90

# See https://github.com/BigDataBiology/SemiBin/issues/155
# A bit of testing showed that scaling tops up at 10~12 processes (on a 16-core
# machine). This leaves a little leeway for a few more processes, but still
# caps it to avoid #155
MAX_PROCESS_ORFS = 24

def reverse_complement(seq: str) -> str:
    '''Reverse complement of a DNA sequence (only uppercase ATGC are complemented)'''
    return seq[::-1].translate(str.maketrans('ATGC', 'TACG'))

codon2aa = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',

    }

def translate(seq: str) -> str:
    '''Translate a DNA sequence into protein (unknown codons become X)'''
    return ''.join([codon2aa.get(seq[i:i+3], 'X') for i in range(0, len(seq), 3)])


def findall(seq: str, pats: list[str]) -> list[int]:
    '''Finds all the matches to the set of patterns given

    Returns
    -------
    positions : list of int
        Sorted start positions of all (possibly overlapping) matches
    '''
    matches = []
    for pat in pats:
        next_start = seq.find(pat)
        while next_start >= 0:
            matches.append(next_start)
            next_start = seq.find(pat, next_start+1)
    matches.sort()
    return matches

def is_atcg(seq: str) -> bool:
    '''Whether the sequence contains only (uppercase) A, T, C, and G'''
    # Empirically this seems to be very fast, even if the code is a bit ridiculous
    return seq.translate(str.maketrans('ATGC', '\n\n\n\n')).count('\n') == len(seq)

def find_orfs_fwd(seq: str, accept_incomplete: bool = False) -> list[ORFPos]:
    '''Find ORFs in the forward strand

    An ORF runs from a start codon to the next in-frame stop codon and must be
    longer than ``MIN_LEN`` and contain only A/C/G/T.

    Parameters
    ----------
    seq : str
        Uppercase DNA sequence
    accept_incomplete : bool, optional
        If True, also report ORFs that are truncated by the start or end of the
        sequence (i.e., missing a start or stop codon)

    Returns
    -------
    orfs : list of ORFPos
    '''
    # orf_start[frame] is the start position of the currently open ORF in
    # that frame, or -1 if there is none
    if accept_incomplete:
        orf_start = [0, 1, 2]
    else:
        orf_start = [-1, -1, -1]
    orfs = []
    positions = findall(seq,
                    START_CODONS + STOP_CODONS)
    for i in positions:
        frame = i % 3
        if seq[i:i+3] in START_CODONS:
            if orf_start[frame] == -1:
                orf_start[frame] = i
        if seq[i:i+3] in STOP_CODONS:
            if orf_start[frame] != -1:
                if i - orf_start[frame] > MIN_LEN:
                    if is_atcg(seq[orf_start[frame]:i+3]):
                        orfs.append(ORFPos(orf_start[frame], i+3, False))
                orf_start[frame] = -1
    if accept_incomplete:
        for frame in range(3):
            if orf_start[frame] != -1:
                if len(seq) - orf_start[frame] > MIN_LEN and is_atcg(seq[orf_start[frame]:]):
                    orfs.append(ORFPos(orf_start[frame], len(seq) - (len(seq)-frame)%3, False))
    return orfs


def find_orfs_rev(seq: str, accept_incomplete: bool = False) -> list[ORFPos]:
    '''Find ORFs in the reverse strand

    Coordinates are given relative to the forward strand (with ``rc=True``).
    See ``find_orfs_fwd`` for the meaning of the arguments.
    '''
    orfs = find_orfs_fwd(reverse_complement(seq), accept_incomplete)
    return [ORFPos(len(seq)-b, len(seq)-a, True) for a,b,_ in orfs][::-1]

def find_orfs(seq: str, accept_incomplete: bool = False) -> list[ORFPos]:
    '''Find ORFs in nucleotide sequence (both strands)

    See ``find_orfs_fwd`` for the meaning of the arguments.
    '''
    return find_orfs_fwd(seq, accept_incomplete) + find_orfs_rev(seq, accept_incomplete)

def extract(seq: str, orf: ORFPos) -> str:
    '''Extract ORF (nucleotide) sequence, reverse complementing if needed'''
    seq = seq[orf.start: orf.end]
    if orf.rc:
        seq = reverse_complement(seq)
    return seq

def orfs_to_fasta(seq: str, header: str, orfs: list[ORFPos]) -> str:
    '''Convert ORFs to (protein) FASTA format

    ORFs are named ``<header>_<index>``, matching prodigal's convention.
    '''
    return ''.join(
            f'>{header}_{i} {orf.start}-{orf.end} {"-1" if orf.rc else "+1"}\n'
            f'{translate(extract(seq, orf))}\n'
                for i,orf in enumerate(orfs))


def get_orfs(header_seq: tuple[str, str]) -> str:
    '''Find ORFs in a single (header, sequence) pair and return them as FASTA

    Takes a tuple so that it can be used with ``Pool.imap``.
    '''
    header, seq = header_seq
    seq = seq.upper()
    orfs = find_orfs(seq, accept_incomplete=True)
    return orfs_to_fasta(seq, header, orfs)


def run_naiveorf(fasta_path: str, num_process: int, output: str) -> str:
    '''Predict ORFs with the built-in naive ORF finder

    Parameters
    ----------
    fasta_path : str
        Input contigs (FASTA)
    num_process : int
        Number of processes to use (capped at ``MAX_PROCESS_ORFS``)
    output : str
        Directory for the output file

    Returns
    -------
    faa_path : str
        Path to the predicted protein sequences (FASTA)
    '''
    oname = path.join(output, 'orfs.faa')
    with open(oname, 'w') as out:
        if num_process == 1:
            for header_seq in fasta.fasta_iter(fasta_path):
                out.write(get_orfs(header_seq))
        else:
            ctx = mp.get_context('spawn')
            with ctx.Pool(processes=min(num_process, MAX_PROCESS_ORFS)) as p:
                outs = p.imap(get_orfs,
                        fasta.fasta_iter(fasta_path),
                        chunksize=8)
                for o in outs:
                    out.write(o)
    return oname
