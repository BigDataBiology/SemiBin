from os import path
import subprocess
import tempfile
import gzip

tmpdir = '/tmp/SemiBin2'
with tempfile.TemporaryDirectory() as tmpdir:
    with open(f'{tmpdir}/S1.fna', 'w') as f:
        f.write('>C1\n')
        for i in range(1000):
            f.write('ATCG')
        f.write('\n')
        f.write('>C2\n')
        for i in range(100):
            f.write('GTAA')
    subprocess.check_call(
        ['SemiBin2', 'split_contigs',
            '-i', f'{tmpdir}/S1.fna',
            '-o', f'{tmpdir}/output',
            '--min-len', '1000',])
    assert path.exists(f'{tmpdir}/output/split_contigs.fna.gz')

    with gzip.open(f'{tmpdir}/output/split_contigs.fna.gz', 'rt') as f:
        assert f.readline() == '>C1_1\n'
        part1 = f.readline()
        assert f.readline() == '>C1_2\n'
        part2 = f.readline()
        assert len(part1) + len(part2) == 4002 # 4000 + 2 newlines
        assert not f.read(1)

