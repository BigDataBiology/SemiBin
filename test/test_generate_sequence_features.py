from SemiBin.main import generate_sequence_features_single, generate_sequence_features_multi
import os
import pytest
import logging
import pandas as pd

def test_generate_seq_feats_multi(tmpdir):
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    os.makedirs(f'{tmpdir}/output_multi',exist_ok=True)
    generate_sequence_features_multi(bams=['test/multi_samples_data/input_multi_sorted1.bam',
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
                         output=f'{tmpdir}/output_multi',
                         contig_fasta='test/multi_samples_data/input_multi.fasta.xz',
                         ratio=0.05,
                         min_length=None,
                         ml_threshold=None,
                         )

    for i in range(10):
        data = pd.read_csv(f'{tmpdir}/output_multi/samples/S{i+1}/data.csv', index_col=0)
        data_split = pd.read_csv(f'{tmpdir}/output_multi/samples/S{i+1}/data_split.csv', index_col=0)
        assert data.shape == (20,146)
        assert data_split.shape == (40,146)



def test_generate_seq_feats_single(tmpdir):
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    os.makedirs(f'{tmpdir}/output_single',exist_ok=True)
    generate_sequence_features_single(
                         bams=['test/single_sample_data/input.sorted.bam'],
                         num_process=1,
                         logger=logger,
                         output=f'{tmpdir}/output_single',
                         contig_fasta='test/single_sample_data/input.fasta',
                         binned_length=2500,
                         must_link_threshold=4000
                         )

    data = pd.read_csv(f'{tmpdir}/output_single/data.csv', index_col=0)
    data_split = pd.read_csv(f'{tmpdir}/output_single/data_split.csv', index_col=0)

    assert data.shape == (40,138)
    assert data_split.shape == (80,136)


def test_generate_seq_feats_coassembly(tmpdir):
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    os.makedirs(f'{tmpdir}/output_coassembly',exist_ok=True)
    generate_sequence_features_single(bams=['test/coassembly_sample_data/input.sorted1.bam',
                               'test/coassembly_sample_data/input.sorted2.bam',
                               'test/coassembly_sample_data/input.sorted3.bam',
                               'test/coassembly_sample_data/input.sorted4.bam',
                               'test/coassembly_sample_data/input.sorted5.bam'],
                         num_process=1,
                         logger=logger,
                         output=f'{tmpdir}/output_coassembly',
                         contig_fasta='test/coassembly_sample_data/input.fasta',
                         binned_length=2500,
                         must_link_threshold=4000
                         )

    data = pd.read_csv(f'{tmpdir}/output_coassembly/data.csv', index_col=0)
    data_split = pd.read_csv(f'{tmpdir}/output_coassembly/data_split.csv', index_col=0)

    assert data.shape == (40,141)
    assert data_split.shape == (80,141)
