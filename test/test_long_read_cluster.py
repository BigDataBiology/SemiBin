'''Tests for the long-read clustering code path'''
import argparse
import logging
from unittest import mock

import numpy as np
import pandas as pd
import pytest
import torch

from SemiBin.long_read_cluster import cluster_long_read
from SemiBin.semi_supervised_model import Semi_encoding_single

N_KMER_FEATURES = 136


class _StopAfterMarkerEstimation(Exception):
    '''Sentinel used to stop cluster_long_read() early'''


# pandas>=3 returns read-only arrays from .values and torch warns about that
@pytest.mark.filterwarnings('ignore:The given NumPy array is not writable')
def test_long_read_fasta_is_closed_before_orffinding(tmp_path):
    '''The concatenated FASTA has to be closed when marker estimation starts

    estimate_seeds() hands the file over to separate processes (the ORF
    finder and HMMER), so anything that is still in the write buffer of this
    process is invisible to them: on most filesystems the input was silently
    truncated and on GPFS-style filesystems it even looked empty.
    '''
    contigs = {'contig1': 'ACGT' * 25, 'contig2': 'TTTT' * 25}
    expected = ''.join(f'>{name}\n{seq}\n' for name, seq in contigs.items())
    seen = {}

    def fake_estimate_seeds(fasta_path, *_args, **_kwargs):
        # this is what prodigal/FragGeneScan/HMMER would get to see
        with open(fasta_path) as f:
            seen['fasta'] = f.read()
        raise _StopAfterMarkerEstimation()

    data = pd.DataFrame(np.zeros((len(contigs), N_KMER_FEATURES + 1)),
                        index=list(contigs))
    args = argparse.Namespace(num_process=1, orf_finder='prodigal',
                              prodigal_output_faa=None)

    with mock.patch('SemiBin.long_read_cluster.estimate_seeds', fake_estimate_seeds):
        with pytest.raises(_StopAfterMarkerEstimation):
            cluster_long_read(logging.getLogger('test'),
                              Semi_encoding_single(N_KMER_FEATURES),
                              data,
                              torch.device('cpu'),
                              False,
                              1,
                              str(tmp_path),
                              contigs,
                              binned_length=0,
                              args=args,
                              minfasta=0)

    assert seen['fasta'] == expected
