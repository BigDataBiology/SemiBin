from collections import namedtuple
import multiprocessing as mp
from os import path
from . import fasta

ORFPos = namedtuple('ORFPos', ['start', 'end', 'rc'])

START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']
MIN_LEN = 90

def reverse_complement(seq):
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

def translate(seq):
    return ''.join([codon2aa.get(seq[i:i+3], 'X') for i in range(0, len(seq), 3)])


def findall(seq, pats):
    '''Finds all the matches to the set of patterns given
    '''
    matches = []
    for pat in pats:
        next_start = seq.find(pat)
        while next_start >= 0:
            matches.append(next_start)
            next_start = seq.find(pat, next_start+1)
    matches.sort()
    return matches

def is_atcg(seq):
    # Empirically this seems to be very fast, even if the code is a bit ridiculous
    return seq.translate(str.maketrans('ATGC', '\n\n\n\n')).count('\n') == len(seq)

def find_orfs_fwd(seq, accept_incomplete=False):
    '''Find ORFs in the forward strand

    accept_incomplete : bool, optional
    '''
    if accept_incomplete:
        active = [0, 1, 2]
    else:
        active = [-1, -1, -1]
    orfs = []
    positions = findall(seq,
                    START_CODONS + STOP_CODONS)
    for i in positions:
        ix = i % 3
        if seq[i:i+3] in START_CODONS:
            if active[ix] == -1:
                active[ix] = i
        if seq[i:i+3] in STOP_CODONS:
            if active[ix] != -1:
                if i - active[ix] > MIN_LEN:
                    if is_atcg(seq[active[ix]:i+3]):
                        orfs.append(ORFPos(active[ix], i+3, False))
                active[ix] = -1
    if accept_incomplete:
        for ix in range(3):
            if active[ix] != -1:
                if len(seq) - active[ix] > MIN_LEN and is_atcg(seq[active[ix]:]):
                    orfs.append(ORFPos(active[ix], len(seq) - (len(seq)-ix)%3, False))
    return orfs


def find_orfs_rev(seq, accept_incomplete=False):
    orfs = find_orfs_fwd(reverse_complement(seq), accept_incomplete)
    return [ORFPos(len(seq)-b, len(seq)-a, True) for a,b,_ in orfs][::-1]

def find_orfs(seq, accept_incomplete=False):
    '''Find ORFs in nucleotide sequence'''
    return find_orfs_fwd(seq, accept_incomplete) + find_orfs_rev(seq, accept_incomplete)

def extract(seq, orf):
    '''Extract ORF sequence'''
    seq = seq[orf.start: orf.end]
    if orf.rc:
        seq = reverse_complement(seq)
    return seq

def orfs_to_fasta(seq, header, orfs):
    '''Convert ORFs to FASTA format'''
    return ''.join(
            f'>{header}_{i} {orf.start}-{orf.end} {"-1" if orf.rc else "+1"}\n'
            f'{translate(extract(seq, orf))}\n'
                for i,orf in enumerate(orfs))
def get_orfs(h_seq):
    h,seq = h_seq
    seq = seq.upper()
    orfs = find_orfs(seq, accept_incomplete=True)
    return orfs_to_fasta(seq, h, orfs)


def run_naiveorf(fasta_path, num_process, output):
    oname = path.join(output, 'orfs.faa')
    with open(oname, 'w') as out:
        if num_process == 1:
            for header,seq in fasta.fasta_iter(fasta_path):
                seq = seq.upper()
                orfs = find_orfs(seq, accept_incomplete=True)
                for i,orf in enumerate(orfs):
                    out.write(f'>{header}_{i} {orf.start}-{orf.end} {"-1" if orf.rc else "+1"}\n')
                    out.write(f'{translate(extract(seq, orf))}\n')
        else:
            ctx = mp.get_context('spawn')
            with ctx.Pool(processes=num_process) as p:
                outs = p.imap(get_orfs,
                        fasta.fasta_iter(fasta_path),
                        chunksize=8)
                for o in outs:
                    out.write(o)
    return oname
