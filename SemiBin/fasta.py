from collections.abc import Iterator
from typing import TextIO


def fasta_iter(
        fname: str | TextIO,
        full_header: bool = False,
        ) -> Iterator[tuple[str, str]]:
    '''Iterate over a (possibly compressed) FASTA file

    Parameters
    ----------
    fname : str or file-like
        Filename or an already open text file. Compression is inferred from
        the filename extension:

        - ``.gz``: gzip
        - ``.bz2``: bzip2
        - ``.xz``: lzma

        Any other extension is read as plain text.
    full_header : bool, optional
        If True, yields the full header line (without the leading ``>``).
        Otherwise (the default), only the first whitespace-delimited word (the
        sequence ID) is returned.

    Yields
    ------
    (header, seq) : tuple of (str, str)
        Sequence lines are concatenated (with whitespace stripped)
    '''
    header = None
    chunks = []
    if hasattr(fname, 'readline'):
        open_fn = lambda f, _ : f
    elif fname.endswith('.gz'):
        import gzip
        open_fn = gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        open_fn = bz2.open
    elif fname.endswith('.xz'):
        import lzma
        open_fn = lzma.open
    else:
        open_fn = open
    with open_fn(fname, 'rt') as f:
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

