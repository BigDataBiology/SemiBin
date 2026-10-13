from __future__ import annotations
import logging
import torch
import pickle
from torch.nn import Linear, LeakyReLU
from torch import nn
from torch.utils.data import Dataset, DataLoader
import os
from .markers import estimate_seeds
from torch.optim import lr_scheduler
import sys


class Semi_encoding_multiple(torch.nn.Module):
    """
    Model for combined features (k-mer and coverage)

    The encoder maps the input to a 100-dimensional embedding. The decoder
    (used only for the reconstruction loss in semi-supervised training) maps
    it back to the input space.

    Parameters
    ----------
    n_features : int
        Number of input features
    """
    def __init__(self, n_features: int):
        super(Semi_encoding_multiple, self).__init__()
        self.encoder1 = torch.nn.Sequential(
            Linear(n_features, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 100),
        )

        self.decoder1 = torch.nn.Sequential(
            Linear(100, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, n_features),
            nn.Sigmoid(),
        )

    def forward(self, input1: torch.Tensor, input2: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        '''Embed a pair of inputs (siamese network)'''
        return self.encoder1(input1), self.encoder1(input2)

    def decoder(self, input1: torch.Tensor, input2: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        '''Reconstruct a pair of inputs from their embeddings'''
        return self.decoder1(input1), self.decoder1(input2)

    def embedding(self, x: torch.Tensor) -> torch.Tensor:
        '''Embed a batch of inputs'''
        return self.encoder1(x)

    def save_with_params_to(self, path: str) -> None:
        '''Save model weights and constructor parameters (see ``model_load``)'''
        torch.save({
                    'model_name': 'Semi_encoding_multiple',
                    'model_state_dict': self.state_dict(),
                    'params': [self.encoder1[0].in_features],
                    }, path)


class Semi_encoding_single(torch.nn.Module):
    """
    Model for k-mer features (coverage is handled outside the model)

    The encoder maps the input to a 100-dimensional embedding. The decoder
    (used only for the reconstruction loss in semi-supervised training) maps
    it back to the input space.

    Parameters
    ----------
    n_features : int
        Number of input features
    """
    def __init__(self, n_features: int):
        super(Semi_encoding_single, self).__init__()
        self.encoder1 = torch.nn.Sequential(
            Linear(n_features, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 100),
        )

        self.decoder1 = torch.nn.Sequential(
            Linear(100, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, n_features),
            nn.Softmax(dim=1),
        )

    def forward(self, input1: torch.Tensor, input2: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        '''Embed a pair of inputs (siamese network)'''
        return self.encoder1(input1), self.encoder1(input2)

    def decoder(self, input1: torch.Tensor, input2: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        '''Reconstruct a pair of inputs from their embeddings'''
        return self.decoder1(input1), self.decoder1(input2)

    def embedding(self, x: torch.Tensor) -> torch.Tensor:
        '''Embed a batch of inputs'''
        return self.encoder1(x)

    def save_with_params_to(self, path: str) -> None:
        '''Save model weights and constructor parameters (see ``model_load``)'''
        torch.save({
                    'model_name': 'Semi_encoding_single',
                    'model_state_dict': self.state_dict(),
                    'params': [self.encoder1[0].in_features],
                    }, path)


def model_load(path: str,
               device: torch.device,
               warn_on_old_format: bool = True,
               ) -> Semi_encoding_single | Semi_encoding_multiple:
    '''Load model from path

    Models saved with ``save_with_params_to`` are loaded safely (with
    ``weights_only=True``). Older models (pickled whole) are loaded with the
    legacy (unsafe) method, logging a warning on how to update them (unless
    ``warn_on_old_format`` is False). Exits the program if loading fails.

    Parameters
    ----------
    path : str
        Model file
    device : torch.device
        Device to load the model to
    warn_on_old_format : bool, optional
        Whether to warn if the model is in the old format

    Returns
    -------
    model : Semi_encoding_single or Semi_encoding_multiple
    '''
    try:
        saved = torch.load(path, map_location=device, weights_only=True)
    except pickle.UnpicklingError:
        import logging
        from time import sleep
        logger = logging.getLogger("SemiBin2")
        if warn_on_old_format:
            logger.warning('There was an error loading the model.')
            logger.warning('This happens if the model was saved by and older version of SemiBin2 as that approach is not compatible with newer versions of PyTorch.')
            logger.warning('We will retry loading the model with the older method, but this may cause errors.')
            new_path = (path.replace('.h5', '.pt') if path.endswith('.h5') else path + '.pt')
            logger.warning('If the model does load, you can resave it in the new format using the following command:\n'
                            f'\n\tSemiBin2 update_model --model {path} --output {new_path}'
                           '\n')
            logger.warning('Then, you can use the new model file in future runs with no warnings (and future-proof).')
        logger.warning('Alright, let\'s try loading the model now using the older approach...')
        sleep(.5)
        try:
            if device == torch.device('cpu'):
                model = torch.load(path, map_location=torch.device('cpu'), weights_only=False)
            else:
                model = torch.load(path, weights_only=False).to(device)
        except Exception as e:
            logger.error(f'Error loading model: {e}')
            logger.error('Your model file is likely incompatible with the current version of PyTorch.')
            logger.error('Please retrain your model or use an older version of PyTorch to convert it.')
            sys.exit(1)
        return model
    if saved['model_name'] == 'Semi_encoding_single':
        model = Semi_encoding_single(saved['params'][0])
    elif saved['model_name'] == 'Semi_encoding_multiple':
        model = Semi_encoding_multiple(saved['params'][0])
    else:
        raise ValueError(f"Unknown model name: '{saved['model_name']}'")
    model.load_state_dict(saved['model_state_dict'])
    return model.to(device)


def loss_function(embedding1: torch.Tensor,
                  embedding2: torch.Tensor,
                  label: torch.Tensor | None,
                  raw_x_1: torch.Tensor,
                  raw_x_2: torch.Tensor,
                  decoder_x_1: torch.Tensor,
                  decoder_x_2: torch.Tensor,
                  is_label: bool = True,
                  ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor] | torch.Tensor:
    '''Loss for semi-supervised training

    Combines a contrastive loss on the embeddings (must-link pairs, with
    ``label == 1``, are pulled together; cannot-link pairs, with
    ``label == 0``, are pushed at least distance 1 apart) and an autoencoder
    reconstruction (MSE) loss.

    Returns
    -------
    If ``is_label`` is True, the tuple ``(loss, supervised_loss,
    unsupervised_loss)``. Otherwise (unlabeled data, ``label`` is ignored),
    only the reconstruction loss.
    '''
    relu = torch.nn.ReLU()
    mse_loss = torch.nn.MSELoss()
    d = torch.norm(embedding1 - embedding2, p=2, dim=1)
    square_pred = torch.square(d)
    margin_square = torch.square(relu(1 - d))

    if is_label:
        supervised_loss = torch.mean(
            label * square_pred + (1 - label) * margin_square)
        unsupervised_loss = 0.5 * \
            mse_loss(decoder_x_1, raw_x_1) + 0.5 * \
            mse_loss(decoder_x_2, raw_x_2)
        loss = supervised_loss + unsupervised_loss
        return loss, supervised_loss, unsupervised_loss

    else:
        unsupervised_loss = 0.5 * \
            mse_loss(decoder_x_1, raw_x_1) + 0.5 * \
            mse_loss(decoder_x_2, raw_x_2)
        return unsupervised_loss


class feature_Dataset(Dataset):
    '''Dataset of labeled pairs of feature vectors

    Item ``i`` is ``(features1[i], features2[i], labels[i])``, where the label
    is 1 for must-link pairs and 0 for cannot-link pairs.
    '''
    def __init__(self, features1, features2, labels):
        self.features1 = features1
        self.features2 = features2
        assert len(features1) == len(features2)
        assert len(features1) == len(labels)
        self.labels = labels

    def __getitem__(self, item):
        return self.features1[item], self.features2[item], self.labels[item]

    def __len__(self) -> int:
        return len(self.features1)


class unsupervised_feature_Dataset(Dataset):
    '''Dataset of unlabeled pairs of feature vectors

    Item ``i`` is ``(features1[i], features2[i])``.
    '''
    def __init__(self, features1, features2):
        self.features1 = features1
        self.features2 = features2
        assert len(features1) == len(features2)

    def __getitem__(self, item):
        return self.features1[item], self.features2[item]

    def __len__(self) -> int:
        return len(self.features1)


def train_semi(logger: logging.Logger,
               out: str,
               contig_fastas: list[str],
               binned_lengths: list[int],
               datas: list[str],
               data_splits: list[str],
               cannot_links: list[str],
               is_combined: bool = True,
               batchsize: int = 2048,
               epochs: int = 20,
               device: torch.device | None = None,
               num_process: int = 8,
               mode: str = 'single',
               orf_finder: str = 'prodigal',
               prodigal_output_faa: str | None = None,
               ) -> Semi_encoding_single | Semi_encoding_multiple:
    """
    Train model from one sample(--mode single) or several samples(--mode several)

    Semi-supervised training uses three types of constraints: must-link
    (halves of the same contig), cannot-link from taxonomic annotation, and
    cannot-link between seed contigs (see ``estimate_seeds``). An autoencoder
    reconstruction loss is also used on all contigs.

    Parameters
    ----------
    logger : logging.Logger
    out : str
        Output directory (marker gene results are saved to
        ``<out>/sample<index>``)
    contig_fastas : list of str
        Contig FASTA file of each sample
    binned_lengths : list of int
        Minimum contig length of each sample
    datas : list of str
        ``data.csv`` file of each sample
    data_splits : list of str
        ``data_split.csv`` file of each sample
    cannot_links : list of str
        Cannot-link file of each sample
    is_combined : bool, optional
        Whether coverage is part of the model input
    batchsize : int, optional
    epochs : int, optional
    device : torch.device, optional
    num_process : int, optional
        Number of threads/processes to use
    mode : str, optional
        'single' or 'several'
    orf_finder : str, optional
        ORF finder to use for marker gene detection
    prodigal_output_faa : str, optional
        Pre-computed ORFs

    Returns
    -------
    model : Semi_encoding_single or Semi_encoding_multiple
        The trained model
    """
    from tqdm import tqdm
    import pandas as pd
    from sklearn.preprocessing import normalize
    import numpy as np
    from .utils import norm_abundance

    train_data = pd.read_csv(datas[0], index_col=0).values
    if not is_combined:
        train_data_input = train_data[:, 0:136]
    else:
        train_data_input = train_data


    torch.set_num_threads(num_process)

    logger.info('Training model...')

    if not is_combined:
        model = Semi_encoding_single(train_data_input.shape[1]).to(device)
    else:
        model = Semi_encoding_multiple(train_data_input.shape[1]).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)

    for epoch in tqdm(range(epochs)):
        for data_index in range(len(contig_fastas)):
            seeds = estimate_seeds(
                contig_fastas[data_index],
                binned_length=binned_lengths[data_index],
                num_process=num_process,
                output=os.path.join(out, 'sample{}'.format(data_index)),
                orf_finder=orf_finder,
                prodigal_output_faa=prodigal_output_faa)
            if epoch == 0:
                logger.debug('Generate training data of {}:'.format(data_index))

            # #generate two inputs
            train_input_1 = []
            train_input_2 = []
            train_labels = []

            # cannot link

            if not os.path.getsize(cannot_links[data_index]):
                logger.error(
                    f"Cannot-link file ({cannot_links[data_index]}) is empty!")
                sys.exit(1)

            cannot_link = pd.read_csv(cannot_links[data_index], sep=',',
                                      header=None).values
            data = pd.read_csv(datas[data_index], index_col=0)
            data.index = data.index.astype(str)
            data_split = pd.read_csv(data_splits[data_index], index_col=0)

            if mode == 'several':
                if data.shape[1] != 138 or data_split.shape[1] != 136:
                    logger.error(
                        "'--train-from-many' training expects single-sample feature files "
                        f"(138 and 136 columns), but got {data.shape[1]} and {data_split.shape[1]} columns.")
                    sys.exit(1)

            contig2ix = {c:ix for ix,c in enumerate(data.index)}
            train_data = data.values
            train_data_must_link = data_split.values

            if not is_combined:
                train_data_input = train_data[:, 0:136]
                train_data_split_input = train_data_must_link
            else:
                if norm_abundance(train_data):
                    norm = np.sum(train_data, axis=0)
                    train_data = train_data / norm
                    train_data_must_link = train_data_must_link / norm
                    train_data_input = normalize(train_data, axis=1, norm='l1')
                    train_data_split_input = normalize(train_data_must_link, axis=1, norm='l1')

                else:
                    train_data_input = train_data
                    train_data_split_input = train_data_must_link

            # cannot link from contig annotation
            for link in cannot_link:
                train_input_1.append(train_data_input[contig2ix[str(link[0])]])
                train_input_2.append(train_data_input[contig2ix[str(link[1])]])
                train_labels.append(0)

            # cannot link from bin seed
            if seeds is not None:
                for i in range(len(seeds)):
                    for j in range(i + 1, len(seeds)):
                        train_input_1.append(
                            train_data_input[contig2ix[str(seeds[i])]])
                        train_input_2.append(
                            train_data_input[contig2ix[str(seeds[j])]])
                        train_labels.append(0)

            # must link from breaking up
            for i in range(0, len(train_data_must_link), 2):
                train_input_1.append(train_data_split_input[i])
                train_input_2.append(train_data_split_input[i + 1])
                train_labels.append(1)

            if epoch == 0:
                logger.debug(
                    'Number of must-link pairs: {}'.format(train_labels.count(1)))
                logger.debug(
                    'Number of cannot-link pairs: {}'.format(train_labels.count(0)))

            dataset = feature_Dataset(train_input_1, train_input_2,
                                      train_labels)
            train_loader = DataLoader(
                dataset=dataset,
                batch_size=batchsize,
                shuffle=True,
                num_workers=0,
                drop_last=True)

            unlabeled_x = train_data_input
            unlabeled_train_input1 = []
            unlabeled_train_input2 = []
            for i in range(len(unlabeled_x)):
                unlabeled_train_input1.append(unlabeled_x[i])
                unlabeled_train_input2.append(unlabeled_x[i])

            dataset_unlabeled = unsupervised_feature_Dataset(
                unlabeled_train_input1, unlabeled_train_input2)
            train_loader_unlabeled = DataLoader(
                dataset=dataset_unlabeled,
                batch_size=batchsize,
                shuffle=True,
                num_workers=0,
                drop_last=True)

            for train_input1, train_input2, train_label in train_loader:
                model.train()
                train_input1 = train_input1.to(device)
                train_input2 = train_input2.to(device)
                train_label = train_label.to(device)
                embedding1, embedding2 = model.forward(
                    train_input1.float(), train_input2.float())
                decoder1, decoder2 = model.decoder(embedding1, embedding2)
                optimizer.zero_grad()
                loss, supervised_loss, unsupervised_loss = loss_function(embedding1.double(), embedding2.double(),
                                                                         train_label.double(), train_input1.double(),
                                                                         train_input2.double(), decoder1.double(),
                                                                         decoder2.double())
                loss = loss.to(device)
                loss.backward()
                optimizer.step()

            for unlabeled_train_input1, unlabeled_train_input2 in train_loader_unlabeled:
                model.train()
                unlabeled_train_input1 = unlabeled_train_input1.to(device)
                unlabeled_train_input2 = unlabeled_train_input2.to(device)
                embedding1, embedding2 = model.forward(
                    unlabeled_train_input1.float(), unlabeled_train_input2.float())
                decoder1, decoder2 = model.decoder(embedding1, embedding2)
                optimizer.zero_grad()
                loss = loss_function(embedding1.double(), embedding2.double(),
                                     None, unlabeled_train_input1.double(),
                                     unlabeled_train_input2.double(), decoder1.double(),
                                     decoder2.double(), is_label=False).to(device)
                loss.backward()
                optimizer.step()

        scheduler.step()

    logger.info('Training finished.')
    return model
