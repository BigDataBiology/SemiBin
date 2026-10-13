from __future__ import annotations
import logging
import torch
from torch.utils.data import DataLoader
from torch.optim import lr_scheduler
import sys
from .semi_supervised_model import Semi_encoding_single, Semi_encoding_multiple, feature_Dataset
from .utils import norm_abundance

def loss_function(embedding1: torch.Tensor, embedding2: torch.Tensor, label: torch.Tensor) -> torch.Tensor:
    '''Contrastive loss

    Must-link pairs (``label == 1``) are pulled together (squared distance),
    while cannot-link pairs (``label == 0``) are pushed at least distance 1
    apart.
    '''
    relu = torch.nn.ReLU()
    d = torch.norm(embedding1 - embedding2, p=2, dim=1)
    square_pred = torch.square(d)
    margin_square = torch.square(relu(1 - d))
    supervised_loss = torch.mean(
        label * square_pred + (1 - label) * margin_square)
    return supervised_loss


def train_self(logger: logging.Logger,
               datapaths: list[str],
               data_splits: list[str],
               is_combined: bool = True,
               batchsize: int = 2048,
               epochs: int = 15,
               device: torch.device | None = None,
               num_process: int = 8,
               mode: str = 'single',
               ) -> Semi_encoding_single | Semi_encoding_multiple:
    """
    Train model from one sample(mode=single) or several samples(mode=several)

    Self-supervised (contrastive) training: must-link pairs are the two halves
    of long contigs (from ``data_split``), while cannot-link pairs are random
    pairs of contigs (resampled every epoch).

    Parameters
    ----------
    logger : logging.Logger
    datapaths : list of str
        ``data.csv`` file of each sample
    data_splits : list of str
        ``data_split.csv`` file of each sample
    is_combined : bool, optional
        Whether coverage is part of the model input
    batchsize : int, optional
    epochs : int, optional
    device : torch.device, optional
    num_process : int, optional
        Number of threads to use
    mode : str, optional
        'single' or 'several'

    Returns
    -------
    model : Semi_encoding_single or Semi_encoding_multiple
        The trained model (it is not saved to disk)
    """
    from tqdm import tqdm
    import pandas as pd
    import numpy as np

    torch.set_num_threads(num_process)

    logger.info('Training model...')

    # Read each (data, data_split) pair once and precompute everything that
    # does not depend on the epoch (CSV parsing, normalization, validation,
    # cannot-link count). The epoch loop below only resamples cannot-link
    # indices and runs the gradient updates.
    prepared = []
    for data_index, (datapath, data_split_path) in enumerate(zip(datapaths, data_splits, strict=True)):
        logger.debug(f'Reading training data for index {data_index}...')
        data = pd.read_csv(datapath, index_col=0)
        data.index = data.index.astype(str)
        data_split = pd.read_csv(data_split_path, index_col=0)

        logger.debug(f'Data shape from file `{datapath}`: {data.shape}')
        logger.debug(f'Data shape from file `{data_split_path}` (split data file): {data_split.shape}')

        if mode == 'several':
            if data.shape[1] != 138 or data_split.shape[1] != 136:
                logger.error(
                    "'--train-from-many' training expects single-sample feature files "
                    f"(138 and 136 columns), but got {data.shape[1]} and {data_split.shape[1]} columns.")
                sys.exit(1)

        train_data = data.values
        train_data_split = data_split.values
        if not is_combined:
            train_data = train_data[:, :136]
        else:
            if norm_abundance(train_data):
                from sklearn.preprocessing import normalize
                norm = np.sum(train_data, axis=0)
                train_data = train_data / norm
                train_data_split = train_data_split / norm
                train_data = normalize(train_data, axis=1, norm='l1')
                train_data_split = normalize(train_data_split, axis=1, norm='l1')

        data_length = len(train_data)
        if data_length == 0:
            logger.error(f'No data for sample {datapath}')
            raise ValueError(f'No data for sample {datapath}')
        elif data_length == 1:
            logger.error(f'Only one data point for sample {datapath} (binning would fail)')
            raise ValueError(f'Only one data point for sample {datapath} (binning would fail)')

        n_must_link = len(train_data_split)
        n_cannot_link = min(n_must_link * 1000 // 2, 4_000_000)
        logger.debug(f'Number of must-link pairs: {n_must_link//2}')
        logger.debug(f'Number of cannot-link pairs: {n_cannot_link}')

        prepared.append((datapath, train_data, train_data_split, data_length, n_cannot_link))

    n_features = prepared[0][1].shape[1]
    if not is_combined:
        model = Semi_encoding_single(n_features)
    else:
        model = Semi_encoding_multiple(n_features)
    model = model.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)

    for epoch in tqdm(range(epochs)):
        for datapath, train_data, train_data_split, data_length, n_cannot_link in prepared:
            # cannot-link data is sampled randomly each epoch
            indices1 = np.random.choice(data_length, size=n_cannot_link)
            indices2 = indices1 + 1 + np.random.choice(data_length - 1,
                                                       size=n_cannot_link)
            indices2 %= data_length

            train_input_1 = np.concatenate(
                                (train_data[indices1],
                                train_data_split[::2]))
            train_input_2 = np.concatenate(
                                    (train_data[indices2],
                                    train_data_split[1::2]))
            train_labels = np.zeros(len(train_input_1), dtype=np.float32)
            train_labels[len(indices1):] = 1
            dataset = feature_Dataset(train_input_1, train_input_2, train_labels)
            train_loader = DataLoader(
                dataset=dataset,
                batch_size=batchsize,
                shuffle=True,
                num_workers=0,
                drop_last=True)

            for train_input1, train_input2, train_label in train_loader:
                model.train()
                train_input1 = train_input1.to(device=device, dtype=torch.float32)
                train_input2 = train_input2.to(device=device, dtype=torch.float32)
                train_label = train_label.to(device=device, dtype=torch.float32)
                embedding1, embedding2 = model.forward(
                    train_input1, train_input2)
                optimizer.zero_grad()
                supervised_loss = loss_function(embedding1.double(), embedding2.double(), train_label.double())
                supervised_loss = supervised_loss.to(device)
                supervised_loss.backward()
                optimizer.step()

        scheduler.step()

    logger.info('Training finished.')
    return model
