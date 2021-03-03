import numpy as np
import os
import subprocess
from atomicwrites import atomic_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys


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

    expect_file(args.contig_fasta)
    expect_file_list(args.bams)
    expect_file_list(args.cannot_link)


def get_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    basepair_sum = 0
    threshold = 0
    whole_len = np.sum(contig_len)
    contig_len.sort(reverse=True)
    index = 0
    while(basepair_sum / whole_len < 0.98):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    threshold = max(threshold, 4000)
    return threshold


def write_bins(namelist, contig_labels, output, contig_dict,
               recluster=False, origin_label=0):
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
        if not recluster:
            if whole_bin_bp >= 200000:
                with atomic_write(os.path.join(output, 'bin.{}.fa'.format(label)), overwrite=True) as ofile:
                    SeqIO.write(bin, ofile, 'fasta')
        else:
            if whole_bin_bp >= 200000:
                with atomic_write(os.path.join(output, 'recluster_{0}.bin.{1}.fa'.format(origin_label, label)), overwrite=True) as ofile:
                    SeqIO.write(bin, ofile, 'fasta')


def cal_kl(m1, m2, v1, v2):
    m1 = max(m1, 1e-6)
    m2 = max(m2, 1e-6)
    v1 = 1 if v1 < 1 else v1
    v2 = 1 if v2 < 1 else v2
    value = np.log(np.sqrt(v2 / v1)) + \
        np.divide(np.add(v1, np.square(m1 - m2)), 2 * v2) - 0.5
    return min(max(value, 1e-6), 1 - 1e-6)


def cal_num_bins(fasta_path, contig_output, hmm_output,
                 seed_output, binned_short):
    if not os.path.exists(contig_output + '.faa'):
        frag_out_log = open(contig_output + '.out', 'w')
        subprocess.check_call(
            ['run_FragGeneScan.pl',
             '-genome={}'.format(fasta_path),
             '-out={}'.format(contig_output),
             '-complete=0',
             '-train=complete',
             '-thread=48',
             ],
            stdout=frag_out_log,
            stderr=subprocess.DEVNULL,
        )

    if not os.path.exists(hmm_output):
        hmm_out_log = open(hmm_output + '.out', 'w')
        subprocess.check_call(
            ['hmmsearch',
             '--domtblout',
             hmm_output,
             '--cut_tc',
             '--cpu', str(48),
             os.path.split(__file__)[0] + '/marker.hmm',
             contig_output + '.faa',
             ],
            stdout=hmm_out_log,
            stderr=subprocess.DEVNULL,
        )

    if not os.path.exists(seed_output):
        if binned_short:
            getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
            subprocess.check_call(
                ['perl', getmarker,
                 hmm_output,
                 fasta_path,
                 str(1001),
                 seed_output,
                 ],
                stdout=open('/dev/null', 'w'),
                stderr=subprocess.DEVNULL,
            )
        else:
            getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
            subprocess.check_call(
                ['perl', getmarker,
                 hmm_output,
                 fasta_path,
                 str(2501),
                 seed_output,
                 ],
                stdout=open('/dev/null', 'w'),
                stderr=subprocess.DEVNULL,
            )
