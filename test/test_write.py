from SemiBin.utils import possibly_compressed_write


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

