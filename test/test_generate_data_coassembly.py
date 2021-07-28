from SemiBin.main import generate_data_single
import os
import pytest
import logging
import pandas as pd

def test_generate_data_coassembly():
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    os.makedirs('output_coassembly',exist_ok=True)
    generate_data_single(bams=['test/coassembly_sample_data/input.sorted1.bam',
                               'test/coassembly_sample_data/input.sorted2.bam',
                               'test/coassembly_sample_data/input.sorted3.bam',
                               'test/coassembly_sample_data/input.sorted4.bam',
                               'test/coassembly_sample_data/input.sorted5.bam'],
                         num_process=1,
                         logger=logger,
                         output='output_coassembly',
                         contig_fasta='test/coassembly_sample_data/input.fasta',
                         binned_length=2500,
                         must_link_threshold=4000
                         )

    data = pd.read_csv('output_coassembly/data.csv',index_col=0)
    data_split = pd.read_csv('output_coassembly/data_split.csv',index_col=0)

    assert data.shape == (40,141)
    assert data_split.shape == (80,141)