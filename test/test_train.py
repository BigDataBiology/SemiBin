from SemiBin.main import training
from SemiBin.fasta import fasta_iter
import os
import pytest
import logging
import pandas as pd

def test_train():
    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    contig_dict = {h:seq for h,seq in fasta_iter('test/train_data/input.fasta')}

    os.makedirs('output_train',exist_ok=True)
    training(contig_fasta = ['test/train_data/input.fasta'],
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
