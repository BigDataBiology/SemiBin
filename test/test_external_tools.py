import os
import stat
import time
import logging
import signal
from unittest import mock

import pytest

from SemiBin.utils import describe_returncode, format_stderr
from SemiBin.orffinding import run_prodigal
from SemiBin.generate_coverage import generate_cov
from SemiBin.markers import estimate_seeds

INPUT_FASTA = 'test/single_sample_data/input.fasta'


def _fake_binary(bindir, name, script):
    os.makedirs(bindir, exist_ok=True)
    path = os.path.join(bindir, name)
    with open(path, 'w') as f:
        f.write('#!/bin/sh\n' + script)
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
    return path


def test_describe_returncode():
    assert describe_returncode(1) == 'exited with return code 1'
    assert describe_returncode(-signal.SIGKILL) == 'was killed by SIGKILL'


def test_format_stderr():
    assert format_stderr('') == ''
    assert format_stderr(b'  \n') == ''
    r = format_stderr(b'first\nsecond\n')
    assert 'first' in r
    assert 'second' in r
    r = format_stderr('\n'.join(f'line{i}' for i in range(100)), max_lines=3)
    assert 'line99' in r
    assert 'line97' in r
    assert 'line96' not in r


def test_run_prodigal_missing(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    with pytest.raises(SystemExit):
        run_prodigal(INPUT_FASTA, 2, str(tmp_path))


def test_run_prodigal_failure_kills_others(tmp_path, monkeypatch, caplog):
    import subprocess
    bindir = str(tmp_path / 'bin')
    # The first chunk fails immediately, the others hang
    _fake_binary(bindir, 'prodigal', '''
case "$2" in
    *contig_0.fa) echo "prodigal: something went wrong" >&2; exit 3 ;;
esac
exec sleep 60
''')
    monkeypatch.setenv('PATH', bindir + os.pathsep + os.environ['PATH'])
    out = tmp_path / 'out'
    out.mkdir()

    processes = []
    real_popen = subprocess.Popen
    def popen(*args, **kwargs):
        p = real_popen(*args, **kwargs)
        processes.append(p)
        return p

    start = time.monotonic()
    with mock.patch('SemiBin.orffinding.subprocess.Popen', side_effect=popen), \
            caplog.at_level(logging.ERROR, logger='SemiBin2'):
        with pytest.raises(SystemExit):
            run_prodigal(INPUT_FASTA, 4, str(out))
    assert time.monotonic() - start < 30
    assert 'exited with return code 3' in caplog.text
    assert 'something went wrong' in caplog.text

    assert len(processes) > 1
    # All processes have been killed & reaped
    for p in processes:
        assert p.returncode is not None


def test_generate_cov_missing_bedtools(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    with pytest.raises(OSError, match='bedtools'):
        generate_cov('input.bam', 0, str(tmp_path), 4000, False, 1000,
                     logging.getLogger('SemiBin2'))


def test_generate_cov_bedtools_failure(tmp_path, monkeypatch):
    bindir = str(tmp_path / 'bin')
    _fake_binary(bindir, 'bedtools', '''
echo "Error: the BAM file is truncated" >&2
exit 1
''')
    monkeypatch.setenv('PATH', bindir + os.pathsep + os.environ['PATH'])
    with pytest.raises(OSError) as e:
        generate_cov('input.bam', 0, str(tmp_path), 4000, False, 1000,
                     logging.getLogger('SemiBin2'))
    assert 'exited with return code 1' in str(e.value)
    assert 'the BAM file is truncated' in str(e.value)
    assert 'input.bam' in str(e.value)


def test_estimate_seeds_missing_hmmsearch(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    with mock.patch('SemiBin.orffinding.run_orffinder') as run_orffinder:
        with pytest.raises(SystemExit):
            estimate_seeds(INPUT_FASTA, 1000, 1, orf_finder='fast-naive')
        assert run_orffinder.call_count == 0


def test_estimate_seeds_hmmsearch_failure(tmp_path, monkeypatch, caplog):
    bindir = str(tmp_path / 'bin')
    _fake_binary(bindir, 'hmmsearch', '''
echo "Error: failed to open the sequence file" >&2
exit 1
''')
    monkeypatch.setenv('PATH', bindir + os.pathsep + os.environ['PATH'])
    with caplog.at_level(logging.ERROR, logger='SemiBin2'):
        with pytest.raises(SystemExit):
            estimate_seeds(INPUT_FASTA, 1000, 1, orf_finder='fast-naive')
    assert 'hmmsearch exited with return code 1' in caplog.text
    assert 'failed to open the sequence file' in caplog.text
    assert INPUT_FASTA in caplog.text
