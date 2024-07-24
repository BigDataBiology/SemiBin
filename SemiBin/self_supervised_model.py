import torch
from torch.utils.data import DataLoader
import os
from torch.optim import lr_scheduler
import sys
from .semi_supervised_model import Semi_encoding_single, Semi_encoding_multiple, feature_Dataset
from .utils import norm_abundance, get_features, normalize_kmer_motif_features

def loss_function(embedding1, embedding2, label):
    relu = torch.nn.ReLU()
    d = torch.norm(embedding1 - embedding2, p=2, dim=1)
    square_pred = torch.square(d)
    margin_square = torch.square(relu(1 - d))
    supervised_loss = torch.mean(
        label * square_pred + (1 - label) * margin_square)
    return supervised_loss




def train_self(logger, out : str, datapaths, data_splits, is_combined=True,
          batchsize=2048, epoches=15, device=None, num_process = 8, mode = 'single'):
    """
    Train model from one sample(mode=single) or several samples(mode=several)

    Saves model to disk and returns it

    Parameters
    ----------
    out : filename to write model to
    """
    from tqdm import tqdm
    import pandas as pd
    from sklearn.preprocessing import normalize
    import numpy as np
    
    train_data = pd.read_csv(datapaths[0], index_col=0)
    features_data = get_features(train_data)
    train_data = train_data.values
    
    features_data_split = get_features(pd.read_csv(data_splits[0], index_col=0))
    
    
    
    if not is_combined:
        train_data = train_data[:, features_data['kmer'] + features_data['motif']]

    torch.set_num_threads(num_process)

    logger.info('Training model...')
    
    if not is_combined:
        model = Semi_encoding_single(train_data.shape[1])
    else:
        model = Semi_encoding_multiple(train_data.shape[1])

    model = model.to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)

    for epoch in tqdm(range(epoches)):
        for data_index, (datapath, data_split_path) in enumerate(zip(datapaths, data_splits)):
            if epoch == 0:
                logger.debug(f'Reading training data for index {data_index}...')

            data = pd.read_csv(datapath, index_col=0)
            data_split = pd.read_csv(data_split_path, index_col=0)

            if mode == 'several':
                if data.shape[1] != (len(features_data['kmer']) + len(features_data['motif']) + 2) or data_split.shape[1] != (len(features_data['kmer']) + len(features_data["motif"])): # + from having a sample column bam_mean + bam_var
                    sys.stderr.write(
                        f"Error: training mode with several only used in single-sample binning!\n")
                    sys.exit(1)

            train_data = data.values
            train_data_split = data_split.values

            if not is_combined:
                train_data = train_data[:, features_data['kmer'] + features_data['motif']]
                if len(features_data["motif"]) > 0:
                    train_data, train_data_split = normalize_kmer_motif_features(train_data, train_data_split)
                
            else:
                if norm_abundance(train_data, features_data):
                    train_data_kmer  = train_data[:, features_data['kmer'] + features_data['motif']]
                    train_data_split_kmer  = train_data_split[:, features_data_split['kmer'] + features_data_split['motif']]
                    
                    if len(features_data["motif"]) > 0:
                        train_data_kmer, train_data_split_kmer = normalize_kmer_motif_features(train_data_kmer, train_data_split_kmer)
                    
                    train_data_depth = train_data[:, features_data['depth']]
                    train_data_depth = normalize(train_data_depth, axis=1, norm='l1')
                    train_data = np.concatenate((train_data_kmer, train_data_depth), axis=1)

                    train_data_split_depth = train_data_split[:, features_data_split['depth']]
                    train_data_split_depth = normalize(train_data_split_depth, axis=1, norm='l1')
                    train_data_split = np.concatenate((train_data_split_kmer, train_data_split_depth), axis = 1)

            
            data_length = len(train_data)
            
            # cannot link data is sampled randomly
            n_cannot_link = min(len(train_data_split) * 1000 // 2, 4_000_000)
            indices1 = np.random.choice(data_length, size=n_cannot_link)
            indices2 = indices1 + 1 + np.random.choice(data_length - 1,
                                                       size=n_cannot_link)
            indices2 %= data_length


            if epoch == 0:
                logger.debug(
                    f'Number of must-link pairs: {len(train_data_split)//2}')
                logger.debug(
                    f'Number of cannot-link pairs: {n_cannot_link}')

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
                # decoder1, decoder2 = model.decoder(embedding1, embedding2)
                optimizer.zero_grad()
                supervised_loss = loss_function(embedding1.double(), embedding2.double(), train_label.double())
                supervised_loss = supervised_loss.to(device)
                supervised_loss.backward()
                optimizer.step()
        scheduler.step()

    logger.info('Training finished.')
    torch.save(model, out)

    return model
