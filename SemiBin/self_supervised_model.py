import torch
from torch.utils.data import DataLoader
import os
from torch.optim import lr_scheduler
import sys
from .semi_supervised_model import Semi_encoding_single, Semi_encoding_multiple, feature_Dataset

def loss_function(embedding1, embedding2, label):
    relu = torch.nn.ReLU()
    d = torch.norm(embedding1 - embedding2, p=2, dim=1)
    square_pred = torch.square(d)
    margin_square = torch.square(relu(1 - d))
    supervised_loss = torch.mean(
        label * square_pred + (1 - label) * margin_square)
    return supervised_loss

def train_self(out, logger, datas, data_splits, is_combined=True,
          batchsize=2048, epoches=15, device=None, num_process = 8, mode = 'single'):
    """
    Train model from one sample(--mode single) or several samples(--mode several)
    """
    from tqdm import tqdm
    import pandas as pd
    from sklearn.preprocessing import normalize
    import numpy as np

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
        for data_index in range(len(datas)):
            if epoch == 0:
                logger.info('Generate training data of {}:'.format(data_index))

            data = pd.read_csv(datas[data_index], index_col=0)
            data.index = data.index.astype(str)
            data_split = pd.read_csv(data_splits[data_index], index_col=0)

            if mode == 'several':
                if data.shape[1] != 138 or data_split.shape[1] != 136:
                    sys.stderr.write(
                        f"Error: training mode with several only used in single-sample binning!\n")
                    sys.exit(1)

            train_data = data.values
            train_data_must_link = data_split.values

            if not is_combined:
                train_data_input = train_data[:, 0:136]
                train_data_split_input = train_data_must_link
            else:
                train_data_input = train_data
                train_data_split_input = train_data_must_link

                if train_data_input.shape[1] - 136 > 20:
                    train_data_kmer = train_data_input[:,0:136]
                    train_data_depth = train_data_input[:, 136:len(data.values[0])]
                    train_data_depth = normalize(train_data_depth, axis=1, norm='l1')
                    train_data_input = np.concatenate((train_data_kmer, train_data_depth), axis=1)

                    train_data_split_kmer = train_data_split_input[:,0:136]
                    train_data_split_depth = train_data_split_input[:, 136:len(data.values[0])]
                    train_data_split_depth = normalize(train_data_split_depth, axis=1, norm='l1')
                    train_data_split_input = np.concatenate((train_data_split_kmer, train_data_split_depth), axis = 1)


            data_length = len(train_data_input)

            n_samples = min(len(train_data_split_input) * 1000 // 2, 4000000)
            indices1 = np.random.choice(data_length, size=n_samples)
            indices2 = indices1 + 1 + np.random.choice(data_length - 1,
                                                       size=n_samples)
            indices2 %= data_length
            train_labels = [0] * n_samples

            train_input_1 = []
            train_input_2 = []

            # must link from breaking up

            for index1, index2 in zip(indices1, indices2):
                train_input_1.append(train_data_input[index1])
                train_input_2.append(train_data_input[index2])

            for i in range(0, len(train_data_must_link), 2):
                train_input_1.append(train_data_split_input[i])
                train_input_2.append(train_data_split_input[i + 1])
                train_labels.append(1)

            if epoch == 0:
                logger.debug(
                    'Number of must-link pairs:{}'.format(train_labels.count(1)))
                logger.debug(
                    'Number of cannot-link pairs:{}'.format(train_labels.count(0)))

            dataset = feature_Dataset(train_input_1, train_input_2,
                                      train_labels)
            train_loader = DataLoader(
                dataset=dataset,
                batch_size=batchsize,
                shuffle=True,
                num_workers=0)

            for train_input1, train_input2, train_label in train_loader:
                model.train()
                train_input1 = train_input1.to(device)
                train_input2 = train_input2.to(device)
                train_label = train_label.to(device)
                embedding1, embedding2 = model.forward(
                    train_input1.float(), train_input2.float())
                # decoder1, decoder2 = model.decoder(embedding1, embedding2)
                optimizer.zero_grad()
                supervised_loss = loss_function(embedding1.double(), embedding2.double(), train_label.double())
                supervised_loss = supervised_loss.to(device)
                supervised_loss.backward()
                optimizer.step()
        scheduler.step()

    logger.info('Training finished.')
    torch.save(model, os.path.join(out, 'model.h5'))

    return model
