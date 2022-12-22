import os
import subprocess
import contextlib
import sys
import shutil
import logging
from .naive_orffinder import run_naiveorf


def run_prodigal(fasta_path, num_process, output):
    from .fasta import fasta_iter

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

    try:
        process = []
        for index in range(next_ix):
            with open(os.path.join(output, f'contig_{index}_log.txt') + '.out', 'w') as prodigal_out_log:
                p = subprocess.Popen(
                        ['prodigal',
                         '-i', os.path.join(output, f'contig_{index}.fa'),
                         '-p', 'meta',
                         '-q',
                         '-m', # See https://github.com/BigDataBiology/SemiBin/issues/87
                         '-a', os.path.join(output, f'contig_{index}.faa')
                         ],
                    stdout=prodigal_out_log)
                process.append(p)

        for p in process:
            p.wait()

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

def run_orffinder(fasta_path, num_process, tdir, orf_finder, prodigal_output_faa):
    '''Run ORF finder (depending on the value or the orf_finder argument'''
    if prodigal_output_faa is not None:
        oname = os.path.join(tdir, 'orfs.faa')
        shutil.copyfile(prodigal_output_faa, oname)
        return oname
    elif orf_finder == 'prodigal':
        return run_prodigal(fasta_path, num_process, tdir)
    elif orf_finder == 'fast-naive':
        logger = logging.getLogger('SemiBin')
        logger.info('Running naive ORF finder')
        return run_naiveorf(fasta_path, num_process, tdir)
    return run_fraggenescan(fasta_path, num_process, tdir)
