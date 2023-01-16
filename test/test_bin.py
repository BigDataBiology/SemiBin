from SemiBin.main import binning, binning_preprocess, binning_long
from SemiBin.cluster import run_embed_infomap, recluster_bins
from SemiBin.fasta import fasta_iter
import os
import logging
import pandas as pd
import numpy as np

def test_bin(tmpdir):
    contig_dict = {h:seq for h,seq in fasta_iter('test/bin_data/input.fasta')}

    odir = f'{tmpdir}/output_test_bin'
    os.makedirs(odir,exist_ok=True)
    binning(num_process=1,
            data='test/bin_data/data.csv',
            max_edges=20,
            max_node=1,
            minfasta=0,
            logger=logging,
            output=odir,
            binned_length=1000,
            device='cpu',
            contig_dict=contig_dict,
            model_path='test/bin_data/model.h5',
            recluster=True,
            random_seed=None,
            environment=None,
            )

    assert len(os.listdir(f'{odir}/output_prerecluster_bins')) > 0
    assert len(os.listdir(f'{odir}/output_recluster_bins')) > 0

    odir = f'{tmpdir}/output_test_bin_long'
    os.makedirs(odir, exist_ok=True)
    binning_long(num_process=1,
            data='test/bin_data/data.csv',
            minfasta=0,
            logger=logging,
            output=odir,
            binned_length=1000,
            device='cpu',
            contig_dict=contig_dict,
            model_path='test/bin_data/model.h5',
            random_seed=None,
            environment=None,
            )

    assert len(os.listdir(f'{odir}/output_bins')) > 0



def test_cluster():
    contig_dict = {h:seq for h,seq in fasta_iter('test/bin_data/input.fasta')}
    is_combined, n_sample, data, model = binning_preprocess(data='test/bin_data/data.csv', depth_metabat2=None, model_path='test/bin_data/model.h5', environment=None, device='cpu')

    _, res = run_embed_infomap(
            logger=logging,
            model=model,
            data=data,
            max_edges=20,
            max_node=1,
            device='cpu',
            is_combined=is_combined,
            num_process=1,
            n_sample=n_sample,
            contig_dict=contig_dict,
            random_seed=None)
    np.testing.assert_array_equal(
            res,
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


def test_recluster():
    contig_dict = {h:seq for h,seq in fasta_iter('test/bin_data/input.fasta')}

    contig_labels = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    reclustered = recluster_bins(logging,
            data=pd.read_csv('test/bin_data/data.csv', index_col=0),
            embedding=np.load('test/bin_data/embedding.npy'),
            n_sample=1,
            is_combined=False,
            contig_labels=contig_labels,
            contig_dict=contig_dict,
            binned_length=1000,
            orf_finder='prodigal',
            num_process=1,
            minfasta=0,
            random_seed=123)

    # Computed with a previous version
    np.testing.assert_array_equal(
            reclustered,
            np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 29, 21, 22, 23, 24, 25, 26,
               26, 28, 27, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 10, 11, 12, 13,
               14, 15, 16, 17, 18, 19]))


