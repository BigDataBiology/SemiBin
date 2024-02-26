import argparse
from atomicwrites import atomic_write
import os

def fasta_iter(fname, full_header=False):
    '''Iterate over a (possibly gzipped) FASTA file

    Parameters
    ----------
    fname : str
        Filename.
            If it ends with .gz, gzip format is assumed
            If .bz2 then bzip2 format is assumed
            if .xz, then lzma format is assumerd
    full_header : boolean (optional)
        If True, yields the full header. Otherwise (the default), only the
        first word

    Yields
    ------
    (h,seq): tuple of (str, str)
    '''
    header = None
    chunks = []
    if hasattr(fname, 'readline'):
        op = lambda f,_ : f
    elif fname.endswith('.gz'):
        import gzip
        op = gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        op = bz2.open
    elif fname.endswith('.xz'):
        import lzma
        op = lzma.open
    else:
        op = open
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                elif full_header:
                    header = line.strip()
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)


def generate_file(contig_file, output, min_length, name):
    os.makedirs(output, exist_ok=True)

    with atomic_write(f"{output}/{name}.fa") as ofile:
        for h, seq in fasta_iter(f"{contig_file}"):
            if len(seq) < min_length:
                continue
            half = len(seq) // 2
            h1 = h + "_1"
            seq1 = seq[:half]
            h2 = h + "_2"
            seq2 = seq[half:]
            ofile.write(f">{h1}\n{seq1}\n")
            ofile.write(f">{h2}\n{seq2}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Generate cannt-link constrains from the output of CAT")
    parser.add_argument('-c', '--contig-file',
                        required=True,
                        help='Path to the contig fasta file corresponding to the CAT output file.',
                        dest='contig_file')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None
                        )
    parser.add_argument('-n', '--name',
                        required=False,
                        help='Name for the splitted contig.',
                        dest='name',
                        default='split'
                        )
    parser.add_argument('-m','--min-length',
                        help='Remove contig whose length smaller than this threshold',
                        required=False,
                        default=0,
                        dest='min_length')


    args = parser.parse_args()
    generate_file(
        args.contig_file,
        args.output,
        args.min_length,
        args.name)


if __name__ == '__main__':
    main()