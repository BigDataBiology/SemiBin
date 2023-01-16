from SemiBin.utils import possibly_compressed_write, write_bins
from SemiBin.fasta import fasta_iter
from glob import glob

def test_write_bins(tmpdir):
    contig_dict = {h:seq for h,seq in fasta_iter('test/bin_data/input.fasta')}
    for comp in ['none', 'gz', 'xz']:
        write_bins(namelist=list(contig_dict.keys()),
            contig_labels=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            contig_seqs=contig_dict,
            output=tmpdir+f'/output_bins_{comp}',
            output_compression=comp)
        ext = ('fa' if comp == 'none' else comp)
        assert len(glob(f'{tmpdir}/output_bins_{comp}/*.{ext}')) == 3


def test_possibly_compressed_write(tmp_path):
    import gzip
    import bz2
    import lzma
    file_content = "Hello World!"
    for compression in [None, "gz", "bz2", "xz"]:
        path = tmp_path / "test.txt"
        if compression is not None:
            path = path.with_suffix(f".{compression}")
        with possibly_compressed_write(str(path)) as f:
            f.write(file_content)
        if compression is None:
            assert path.read_text() == file_content
        elif compression == "gz":
            assert gzip.open(path, 'rt').read() == file_content
        elif compression == "bz2":
            assert bz2.open(path, 'rt').read() == file_content
        elif compression == "xz":
            assert lzma.open(path, 'rt').read() == file_content
        else:
            raise ValueError(f"Unknown compression: {compression}")

