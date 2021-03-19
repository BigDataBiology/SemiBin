from s3n2bin.main import binning
import os
import pytest
import logging
import pandas as pd
from Bio import SeqIO

def test_bin():
    logger = logging.getLogger('S3N2Bin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    contig_length_dict = {}
    contig_dict = {}
    handle = 'test/bin_data/input.fasta'
    for seq_record in SeqIO.parse(handle, "fasta"):
        contig_length_dict[str(seq_record.id).strip(
            '')] = len((seq_record.seq))
        contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

    os.makedirs('output_bin',exist_ok=True)
    binning(contig_fasta='test/bin_data/input.fasta',
            bams=['test/bin_data/input.sorted.bam'],
            num_process=1,
            data='test/bin_data/data.csv',
            data_split='test/bin_data/data_split.csv',
            cannot_link='test/bin_data/cannot.txt',
            batchsize=2048,
            epoches=1,
            max_edges=20,
            max_node=1,
            minfasta=0,
            logger=logger,
            output='output_bin',
            binned_short=True,
            device='cpu',
            contig_length_dict=contig_length_dict,
            contig_dict=contig_dict
            )

    assert len(os.listdir('output_bin/output_bins')) > 0
    assert len(os.listdir('output_bin/output_recluster_bins')) > 0