from SemiBin.main import generate_data_single
import os
import pytest
import logging
import pandas as pd

def test_generate_data_single():
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    os.makedirs('output_single',exist_ok=True)
    generate_data_single(bams=['test/single_sample_data/input.sorted.bam'],
                         num_process=1,
                         logger=logger,
                         output='output_single',
                         contig_fasta='test/single_sample_data/input.fasta',
                         binned_length=2500,
                         must_link_threshold=4000
                         )

    data = pd.read_csv('output_single/data.csv',index_col=0)
    data_split = pd.read_csv('output_single/data_split.csv',index_col=0)

    assert data.shape == (40,138)
    assert data_split.shape == (80,136)



