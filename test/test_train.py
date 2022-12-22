from SemiBin.main import training
from SemiBin.fasta import fasta_iter
import os
import logging
import pandas as pd

def test_train(tmpdir):
    contig_dict = {h:seq for h,seq in fasta_iter('test/train_data/input.fasta')}

    odir = f'{tmpdir}/output_train'
    os.makedirs(odir)
    training(contig_fasta = ['test/train_data/input.fasta'],
            num_process = 1,
            data = ['test/train_data/data.csv'],
            data_split = ['test/train_data/data_split.csv'],
            cannot_link = ['test/train_data/cannot.txt'],
            batchsize = 2048,
            epoches = 1,
            logger = logging,
            output = odir,
            device = 'cpu',
            mode = 'single',
            ratio=0.05,
            min_length=None,
            training_mode = 'semi'
            )

    assert os.path.exists(f'{odir}/model.h5')

def test_train_self(tmpdir):
    contig_dict = {h:seq for h,seq in fasta_iter('test/train_data/input.fasta')}
    odir = f'{tmpdir}/output_train_self'
    os.makedirs(odir)
    training(contig_fasta = ['test/train_data/input.fasta'],
            num_process = 1,
            data = ['test/train_data/data.csv'],
            data_split = ['test/train_data/data_split.csv'],
            cannot_link = ['test/train_data/cannot.txt'],
            batchsize = 2048,
            epoches = 1,
            logger = logging,
            output = odir,
            device = 'cpu',
            mode = 'single',
            ratio=0.05,
            min_length=None,
            training_mode = 'self'
            )

    assert os.path.exists(f'{odir}/model.h5')
