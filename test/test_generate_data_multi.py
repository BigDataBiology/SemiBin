from SemiBin.main import generate_data_multi
import os
import pytest
import logging
import pandas as pd

def test_generate_data_multi():
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    os.makedirs('output_multi',exist_ok=True)
    generate_data_multi(bams=['test/multi_samples_data/input_multi_sorted1.bam',
                              'test/multi_samples_data/input_multi_sorted2.bam',
                              'test/multi_samples_data/input_multi_sorted3.bam',
                              'test/multi_samples_data/input_multi_sorted4.bam',
                              'test/multi_samples_data/input_multi_sorted5.bam',
                              'test/multi_samples_data/input_multi_sorted6.bam',
                              'test/multi_samples_data/input_multi_sorted7.bam',
                              'test/multi_samples_data/input_multi_sorted8.bam',
                              'test/multi_samples_data/input_multi_sorted9.bam',
                              'test/multi_samples_data/input_multi_sorted10.bam'],
                         num_process=1,
                         separator=':',
                         logger=logger,
                         output='output_multi',
                         contig_fasta='test/multi_samples_data/input_multi.fasta',
                         ratio=0.05,
                         min_length=None,
                         ml_threshold=None,
                         )

    for i in range(10):
        data = pd.read_csv('output_multi/samples/S{}/data.csv'.format(i+1),index_col=0)
        data_split = pd.read_csv('output_multi/samples/S{}/data_split.csv'.format(i+1),index_col=0)
        assert data.shape == (20,146)
        assert data_split.shape == (40,146)
