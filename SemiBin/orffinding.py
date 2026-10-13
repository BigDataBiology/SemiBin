import os
import subprocess
import contextlib
import sys
import shutil
import logging

from .utils import describe_returncode, format_stderr


def run_prodigal(fasta_path: str, num_process: int, output: str) -> str:
    '''Predict ORFs with prodigal (in metagenome mode)

    The input is split into ``num_process`` chunks of approximately equal size
    (in basepairs) and a prodigal process is run on each chunk in parallel.

    Exits the program (``sys.exit(1)``) if prodigal is not available or fails.

    Parameters
    ----------
    fasta_path : str
        Input contigs (FASTA)
    num_process : int
        Number of prodigal processes to run in parallel
    output : str
        Directory for intermediate and output files

    Returns
    -------
    faa_path : str
        Path to the predicted protein sequences (FASTA)
    '''
    from .fasta import fasta_iter
    logger = logging.getLogger('SemiBin2')
    logger.info('Running prodigal')

    contigs = {}
    for h, seq in fasta_iter(fasta_path):
        contigs[h] = seq

    total_len = sum(len(s) for s in contigs.values())
    chunk_target_len = total_len // num_process

    cur_chunk_len = chunk_target_len + 1
    n_chunks = 0
    out = None
    with contextlib.ExitStack() as stack:
        for h, seq in contigs.items():
            if cur_chunk_len > chunk_target_len and n_chunks < num_process:
                if out is not None:
                    out.close()
                out = open(os.path.join(output, 'contig_{}.fa'.format(n_chunks)), 'wt')
                out = stack.enter_context(out)

                cur_chunk_len = 0
                n_chunks += 1
            out.write(f'>{h}\n{seq}\n')
            cur_chunk_len += len(seq)

    prodigal_path = shutil.which('prodigal')
    if prodigal_path is None:
        logger.error('prodigal binary not found in PATH. Please install prodigal '
                     'or use a different ORF finder (`--orf-finder=fast-naive` or `--orf-finder=fraggenescan`)')
        sys.exit(1)

    processes = []
    try:
        logger.debug(f'Running prodigal on {n_chunks} FASTA files')
        for index in range(n_chunks):
            with open(os.path.join(output, f'contig_{index}_log.txt') + '.out', 'w') as prodigal_out_log, \
                    open(os.path.join(output, f'contig_{index}_log.txt') + '.err', 'w') as prodigal_err_log:
                p = subprocess.Popen(
                        [prodigal_path,
                         '-i', os.path.join(output, f'contig_{index}.fa'),
                         '-p', 'meta',
                         '-q',
                         '-m', # See https://github.com/BigDataBiology/SemiBin/issues/87
                         '-a', os.path.join(output, f'contig_{index}.faa')
                         ],
                    stdout=prodigal_out_log,
                    stderr=prodigal_err_log)
                processes.append(p)

        for index, p in enumerate(processes):
            rc = p.wait()
            if rc != 0:
                with open(os.path.join(output, f'contig_{index}_log.txt') + '.err') as prodigal_err_log:
                    err = format_stderr(prodigal_err_log.read())
                logger.error(f'prodigal {describe_returncode(rc)} '
                             f'(when processing chunk {index + 1} of {n_chunks} of `{fasta_path}`).{err}')
                sys.exit(1)

    except OSError as e:
        logger.error(f"Running prodigal failed: {e}")
        sys.exit(1)

    else:
        logger.debug('All prodigal processes finished successfully')

    finally:
        # Ensure that no prodigal process is left running if we are exiting
        # early (because one of them failed or we were interrupted)
        for p in processes:
            if p.poll() is None:
                p.kill()
        for p in processes:
            p.wait()

    contig_output = os.path.join(output, 'contigs.faa')
    with open(contig_output, 'w') as f:
        for index in range(n_chunks):
            with open(os.path.join(output, 'contig_{}.faa'.format(index)), 'r') as infile:
                f.write(infile.read())
    return contig_output


def run_fraggenescan(fasta_path: str, num_process: int, output: str) -> str:
    '''Predict ORFs with FragGeneScan

    Exits the program (``sys.exit(1)``) if FragGeneScan is not available or
    fails.

    Parameters
    ----------
    fasta_path : str
        Input contigs (FASTA)
    num_process : int
        Number of threads
    output : str
        Directory for output files

    Returns
    -------
    faa_path : str
        Path to the predicted protein sequences (FASTA)
    '''
    logger = logging.getLogger('SemiBin2')
    try:
        contig_output = os.path.join(output, 'contigs.faa')
        with open(contig_output + '.out', 'w') as frag_out_log:
            # We need to call FragGeneScan instead of the Perl wrapper because the
            # Perl wrapper does not handle filepaths correctly if they contain spaces
            # This binary does not handle return codes correctly, though, so we
            # cannot use `check_call`:
            fgs_path = shutil.which('FragGeneScan')
            if fgs_path is None:
                raise SystemError("FragGeneScan binary not found in PATH")
            subprocess.call(
                [fgs_path,
                 '-s', fasta_path,
                 '-o', contig_output,
                 '-w', str(0),
                 '-t', 'complete',
                 '-p', str(num_process),
                 ],
                stdout=frag_out_log,
            )
    except Exception as e:
        logger.error(f"Running FragGeneScan failed: {e}")
        sys.exit(1)
    return contig_output + '.faa'


def run_orffinder(
        fasta_path: str,
        num_process: int,
        tdir: str,
        orf_finder: str,
        prodigal_output_faa: str | None) -> str:
    '''Run the ORF finder selected by the ``orf_finder`` argument

    Parameters
    ----------
    fasta_path : str
        Input contigs (FASTA)
    num_process : int
        Number of processes to use
    tdir : str
        Directory for intermediate and output files
    orf_finder : str
        One of 'prodigal', 'fast-naive', or 'fraggenescan'
    prodigal_output_faa : str or None
        If not None, pre-computed ORFs (protein FASTA) which are used instead
        of running an ORF finder

    Returns
    -------
    faa_path : str
        Path to the predicted protein sequences (FASTA)
    '''
    if prodigal_output_faa is not None:
        oname = os.path.join(tdir, 'orfs.faa')
        shutil.copyfile(prodigal_output_faa, oname)
        return oname
    elif orf_finder == 'prodigal':
        return run_prodigal(fasta_path, num_process, tdir)
    elif orf_finder == 'fast-naive':
        from .naive_orffinder import run_naiveorf
        logger = logging.getLogger('SemiBin2')
        logger.info('Running naive ORF finder')
        return run_naiveorf(fasta_path, num_process, tdir)
    elif orf_finder == 'fraggenescan':
        return run_fraggenescan(fasta_path, num_process, tdir)
    else:
        logger = logging.getLogger('SemiBin2')
        logger.error(f"Unknown ORF finder {orf_finder}")
        sys.exit(1)
