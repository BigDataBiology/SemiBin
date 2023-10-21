import os
import subprocess
import tempfile
import gzip

with tempfile.TemporaryDirectory() as tmpdir:
    with open(f'{tmpdir}/S1.fna', 'w') as f:
        f.write('>C1\n')
        f.write('ATCG\n')
        f.write('>C2\n')
        f.write('AAAA\n')
    with open(f'{tmpdir}/S2.fna', 'w') as f:
        f.write('>C1\n')
        f.write('CGAT\n')
        f.write('>C2\n')
        f.write('TTTT\n')
    subprocess.check_call(
        ['SemiBin2', 'concatenate_fasta',
            '-i', f'{tmpdir}/S1.fna', f'{tmpdir}/S2.fna',
            '-o', f'{tmpdir}/output'])
    assert os.path.exists(f'{tmpdir}/output/concatenated.fa.gz')
    with gzip.open(f'{tmpdir}/output/concatenated.fa.gz', 'rt') as f:
        assert f.readline() == '>S1:C1\n'
        assert f.readline() == 'ATCG\n'
        assert f.readline() == '>S1:C2\n'
        assert f.readline() == 'AAAA\n'
        assert f.readline() == '>S2:C1\n'
        assert f.readline() == 'CGAT\n'
        assert f.readline() == '>S2:C2\n'
        assert f.readline() == 'TTTT\n'
        assert not f.read(1)

