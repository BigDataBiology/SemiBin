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


# https://github.com/BigDataBiology/SemiBin/issues/137
def test_regression_137_semi(tmpdir):
    from SemiBin.semi_supervised_model import train
    odir = f'{tmpdir}/output_train_semi'
    os.makedirs(odir)
    # 40 elements plus header: 41
    assert len(open('test/train_data/data.csv', 'r').readlines()) == 41
    model = train(odir,
                       contig_fastas = ['test/train_data/input.fasta'],
                       logger = logging,
                       binned_lengths = [1000],
                       datas=['test/train_data/data.csv'],
                       data_splits=['test/train_data/data_split.csv'],
                       cannot_links = ['test/train_data/cannot.txt'],
                       is_combined=False,
                       batchsize=39, # test/semi_data/data.csv has 40 elements so bug is triggered when batchsize is 39
                       epoches=1,
                       device='cpu',
                       num_process=1,
                       mode='single',
                       )
    assert os.path.exists(f'{odir}/model.h5')


# https://github.com/BigDataBiology/SemiBin/issues/137
def test_regression_137_self(tmpdir):
    from SemiBin.self_supervised_model import train_self
    odir = f'{tmpdir}/output_train_self'
    os.makedirs(odir)
    # 40 elements plus header: 41
    assert len(open('test/train_data/data.csv', 'r').readlines()) == 41
    # 80 elements plus header: 81
    assert len(open('test/train_data/data_split.csv', 'r').readlines()) == 81

    # Training adds len(<split>) * 1000//2 + 40 so that the total data is 40040
    # To trigger the bug, batchsize is set to 40039
    model = train_self(out = odir,
                       logger = logging,
                       datapaths=['test/train_data/data.csv'],
                       data_splits=['test/train_data/data_split.csv'],
                       is_combined=False,
                       batchsize=40039,
                       epoches=1,
                       device='cpu',
                       num_process=1,
                       mode='single',
                       )
    assert os.path.exists(f'{odir}/model.h5')

