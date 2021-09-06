from SemiBin.main import binning
from SemiBin.fasta import fasta_iter
import os
import pytest
import logging
import pandas as pd

def test_bin():
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    contig_dict = {h:seq for h,seq in fasta_iter('test/bin_data/input.fasta')}

    os.makedirs('output_bin',exist_ok=True)
    binning(num_process=1,
            data='test/bin_data/data.csv',
            max_edges=20,
            max_node=1,
            minfasta=0,
            logger=logger,
            output='output_bin',
            binned_length=1000,
            device='cpu',
            contig_dict=contig_dict,
            model_path='test/bin_data/model.h5',
            recluster=True,
            random_seed=None,
            environment=None,
            )

    assert len(os.listdir('output_bin/output_bins')) > 0
    assert len(os.listdir('output_bin/output_recluster_bins')) > 0
