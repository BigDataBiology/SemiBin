import tempfile
from .utils import maybe_uncompress, fasta_iter
import os
import sys
import subprocess

__normalize_marker_trans = {
    'TIGR00388': 'TIGR00389',
    'TIGR00471': 'TIGR00472',
    'TIGR00408': 'TIGR00409',
    'TIGR02386': 'TIGR02387',
}

def get_marker(hmmout,
               fasta_path = None,
               min_contig_len = None,
               multi_mode: bool = False,
               orf_finder = None,
               contig_to_marker = False):
    '''Parse HMM output file and return markers
    '''
    import pandas as pd
    data = pd.read_table(hmmout, sep=r'\s+',  comment='#', header=None,
                         usecols=(0,3,5,15,16), names=['orf', 'gene', 'qlen', 'qstart', 'qend'])
    if not len(data):
        return []
    data['gene'] = data['gene'].map(lambda m: __normalize_marker_trans.get(m , m))
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

def estimate_seeds(fasta_path,
                 binned_length,
                 num_process,
                 multi_mode: bool = False,
                 output = None,
                 orf_finder: str = 'prodigal',
                 prodigal_output_faa = None):
    '''Estimate number of bins from a FASTA file

    Parameters
    fasta_path: path
    binned_length: int (minimal contig length)
    num_process: int (number of CPUs to use)
    multi_mode: bool, optional (if True, treat input as resulting from concatenating multiple files)

    Returns
    -------
    '''
    from .orffinding import run_orffinder
    with tempfile.TemporaryDirectory() as tdir:
        fasta_path = maybe_uncompress(fasta_path, tdir)
        if output is not None:
            if os.path.exists(os.path.join(output, 'markers.hmmout')):
                return get_marker(os.path.join(output, 'markers.hmmout'),
                                    fasta_path, binned_length, multi_mode, orf_finder=orf_finder)
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
        except Exception as e:
            if os.path.exists(hmm_output):
                os.remove(hmm_output)
            sys.stderr.write(
                    f"Error: Running hmmsearch failed: {e}\n")
            sys.exit(1)

        return get_marker(hmm_output, fasta_path, binned_length, multi_mode, orf_finder=orf_finder)

