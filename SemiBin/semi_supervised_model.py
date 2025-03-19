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
    Model for combined features
    """
    def __init__(self, num):
        super(Semi_encoding_multiple, self).__init__()
        self.encoder1 = torch.nn.Sequential(
            Linear(num, 512),
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
            Linear(512, num),
            nn.Sigmoid(),
        )

    def forward(self, input1, input2):
        return self.encoder1(input1), self.encoder1(input2)

    def decoder(self, input1, input2):
        return self.decoder1(input1), self.decoder1(input2)

    def embedding(self, input):
        return self.encoder1(input)

    def save_with_params_to(self, path):
        torch.save({
                    'model_name': 'Semi_encoding_multiple',
                    'model_state_dict': self.state_dict(),
                    'params': [self.encoder1[0].in_features],
                    }, path)


class Semi_encoding_single(torch.nn.Module):
    """
    Model for k-mer features
    """
    def __init__(self, num):
        super(Semi_encoding_single, self).__init__()
        self.encoder1 = torch.nn.Sequential(
            Linear(num, 512),
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
            Linear(512, num),
            nn.Softmax(dim=1),
        )

    def forward(self, input1, input2):
        return self.encoder1(input1), self.encoder1(input2)

    def decoder(self, input1, input2):
        return self.decoder1(input1), self.decoder1(input2)

    def embedding(self, input):
        return self.encoder1(input)

    def save_with_params_to(self, path):
        torch.save({
                    'model_name': 'Semi_encoding_single',
                    'model_state_dict': self.state_dict(),
                    'params': [self.encoder1[0].in_features],
                    }, path)


def model_load(path, device, warn_on_old_format=True):
    '''Load model from path'''
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
                model = torch.load(path, map_location=torch.device('cpu'))
            else:
                model = torch.load(path).to(device)
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
    model.load_state_dict(saved['model_state_dict'])
    return model.to(device)


def loss_function(embedding1, embedding2, label, raw_x_1,
                  raw_x_2, decoder_x_1, decoder_x_2, is_label=True):
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
    def __init__(self, embedding1, embedding2, labels):
        self.embedding1 = embedding1
        self.embedding2 = embedding2
        assert len(embedding1) == len(embedding2)
        assert len(embedding1) == len(labels)
        self.labels = labels

    def __getitem__(self, item):
        return self.embedding1[item], self.embedding2[item], self.labels[item]

    def __len__(self):
        return len(self.embedding1)


class unsupervised_feature_Dataset(Dataset):
    def __init__(self, embedding1, embedding2):
        self.embedding1 = embedding1
        self.embedding2 = embedding2
        assert len(embedding1) == len(embedding2)

    def __getitem__(self, item):
        return self.embedding1[item], self.embedding2[item]

    def __len__(self):
        return len(self.embedding1)


def train_semi(logger, out, contig_fastas, binned_lengths, datas, data_splits, cannot_links, is_combined=True,
          batchsize=2048, epoches=20, device=None, num_process = 8, mode = 'single', orf_finder = 'prodigal',
          prodigal_output_faa=None):
    """
    Train model from one sample(--mode single) or several samples(--mode several)
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

    for epoch in tqdm(range(epoches)):
        for data_index in range(len(contig_fastas)):
            seed = estimate_seeds(
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
                sys.stderr.write(
                    f"Error: Cannot-link file ({cannot_links[data_index]}) is empty!\n")
                sys.exit(1)

            cannot_link = pd.read_csv(cannot_links[data_index], sep=',',
                                      header=None).values
            data = pd.read_csv(datas[data_index], index_col=0)
            data.index = data.index.astype(str)
            data_split = pd.read_csv(data_splits[data_index], index_col=0)

            if mode == 'several':
                if data.shape[1] != 138 or data_split.shape[1] != 136:
                    sys.stderr.write(
                        f"Error: training mode with several only used in single-sample binning!\n")
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
            if seed is not None:
                for i in range(len(seed)):
                    for j in range(i + 1, len(seed)):
                        train_input_1.append(
                            train_data_input[contig2ix[str(seed[i])]])
                        train_input_2.append(
                            train_data_input[contig2ix[str(seed[j])]])
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
                    unlabeled_train_input1.float(), unlabeled_train_input1.float())
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
