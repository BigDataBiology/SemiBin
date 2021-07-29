from SemiBin.main import training
import os
import pytest
import logging
import pandas as pd
from Bio import SeqIO

def test_train():
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    contig_length_dict = {}
    contig_dict = {}
    handle = 'test/train_data/input.fasta'
    for seq_record in SeqIO.parse(handle, "fasta"):
        contig_length_dict[str(seq_record.id).strip(
            '')] = len((seq_record.seq))
        contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

    os.makedirs('output_train',exist_ok=True)
    training(contig_fasta = ['test/train_data/input.fasta'],
            bams = ['test/train_data/input.sorted.bam'],
            num_process = 1,
            data = ['test/train_data/data.csv'],
            data_split = ['test/train_data/data_split.csv'],
            cannot_link = ['test/train_data/cannot.txt'],
            batchsize = 2048,
            epoches = 1,
            logger = logger,
            output = 'output_train',
            device = 'cpu',
            mode = 'single',
            ratio=0.05,
            min_length=None,
            )

    assert os.path.exists('output_train/model.h5')