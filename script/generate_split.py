import argparse
from atomicwrites import atomic_write
from SemiBin.fasta import fasta_iter
import os

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