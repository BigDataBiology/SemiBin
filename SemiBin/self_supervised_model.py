import torch
from torch.utils.data import DataLoader
import os
from torch.optim import lr_scheduler
import sys
from .semi_supervised_model import Semi_encoding_single, Semi_encoding_multiple, feature_Dataset
from .utils import norm_abundance

def loss_function(embedding1, embedding2, label):
    relu = torch.nn.ReLU()
    d = torch.norm(embedding1 - embedding2, p=2, dim=1)
    square_pred = torch.square(d)
    margin_square = torch.square(relu(1 - d))
    supervised_loss = torch.mean(
        label * square_pred + (1 - label) * margin_square)
    return supervised_loss


def check_motif(column):
    """
    Check if a column is a motif.
    
    Parameters:
        column (str): The column name to check.
    
    Returns:
        bool: True if the column is a motif, False otherwise.
    """
    try:
        motif, mod_pos = column.split('_')
        mod, pos = mod_pos.split('-')
        if mod in ["m", "a"] and int(pos) in range(0, 15):
            return True
    except:
        return False
    
    
    return column.startswith('motif_')

def get_features(df):
    """
    Takes a DataFrame and extracts indices of features to populate the provided features dictionary.
    Specific features are extracted based on the column names:
    - Indices of columns ending in 'bam_mean' or 'bam_var' are considered depth features.
    
    Parameters:
        df (pd.DataFrame): The DataFrame from which to extract features.
        features_dict (dict): The dictionary to populate with feature indices.
    """
    features_dict = {
        'kmer': list(range(136)),
        'depth': [],
        'motif': []
    }
    
    try:
        columns = df.columns
        # Populate 'depth' with indices of columns ending with 'bam_mean' or 'bam_var'
        features_dict['depth'] = [i for i, column in enumerate(columns) if column.endswith('bam_mean') or column.endswith('bam_var')]
        
        # Populate 'motif' with indices of columns
        features_dict['motif'] = [i for i, column in enumerate(columns) if check_motif(column)]
        
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return features_dict


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
                if data.shape[1] != len(features_data['kmer'] + features_data['motif'] + 2) or data_split.shape[1] != len(features_data['kmer'] + features_data["motif"]): # + from having a sample column bam_mean + bam_var
                    sys.stderr.write(
                        f"Error: training mode with several only used in single-sample binning!\n")
                    sys.exit(1)

            train_data = data.values
            train_data_split = data_split.values

            if not is_combined:
                train_data = train_data[:, features_data['kmer'] + features_data['motif']]
            else:
                if norm_abundance(train_data, features_data):
                    train_data_kmer  = train_data[:, features_data['kmer'] + features_data['motif']]
                    train_data_depth = train_data[:, features_data['depth']]
                    train_data_depth = normalize(train_data_depth, axis=1, norm='l1')
                    train_data = np.concatenate((train_data_kmer, train_data_depth), axis=1)

                    train_data_split_kmer  = train_data_split[:, features_data_split['kmer'] + features_data_split['motif']]
                    train_data_split_depth = train_data_split[:, features_data_split['depth']]
                    train_data_split_depth = normalize(train_data_split_depth, axis=1, norm='l1')
                    train_data_split = np.concatenate((train_data_split_kmer, train_data_split_depth), axis = 1)

            
            # informed cannot-links
            ## Calculate kmer cosine distance
            from sklearn.metrics.pairwise import cosine_distances
            kmer_profiles = train_data[:, features_data['kmer']]
            kmer_distances = cosine_distances(kmer_profiles)

            ## Calculate motif hamming distance
            methylation_profiles = train_data[:, features_data['motif']]
            methylation_binary = (methylation_profiles > 0.5).astype(int)
            
            ### Calculate the Hamming distance using broadcasting
            motif_distances = np.sum(methylation_binary[:, np.newaxis] != methylation_binary[np.newaxis, :], axis=2)

            ## Combine values
            from sklearn.preprocessing import MinMaxScaler
            # Normalize the matrices
            cosine_distances_normalized = MinMaxScaler().fit_transform(kmer_distances)
            hamming_distances_normalized = MinMaxScaler().fit_transform(motif_distances)

            # Define weights for the distance matrices
            w1 = 0.5  # Weight for cosine distance
            w2 = 0.5  # Weight for Hamming distance

            # Combine the distance matrices
            combined_distance_matrix = w1 * cosine_distances_normalized + w2 * hamming_distances_normalized

            # Get the upper triangle indices of the combined distance matrix to avoid duplicate pairs
            upper_triangle_indices = np.triu_indices_from(combined_distance_matrix, k=1)

            # Get the combined distances for these pairs
            distances = combined_distance_matrix[upper_triangle_indices]

            # Filter pairs where the distance is at least 0.1
            valid_pairs = np.where(distances >= 0.1)[0]

            # Get the corresponding indices
            filtered_indices1 = upper_triangle_indices[0][valid_pairs]
            filtered_indices2 = upper_triangle_indices[1][valid_pairs]
            # Calculate the number of cannot-links
            n_cannot_link = min(len(data_split) * 1000 // 2, 4_000_000)

            # Ensure we do not exceed the available number of valid pairs
            n_cannot_link = min(n_cannot_link, len(filtered_indices1))

            # Randomly sample indices from the filtered pairs
            sample_indices = np.random.choice(len(filtered_indices1), size=n_cannot_link, replace=False)

            # Get the sampled cannot-link indices
            indices1 = filtered_indices1[sample_indices]
            indices2 = filtered_indices2[sample_indices]
            
            data_length = len(train_data)
            
            # cannot link data is sampled randomly
            # n_cannot_link = min(len(train_data_split) * 1000 // 2, 4_000_000)
            # indices1 = np.random.choice(data_length, size=n_cannot_link)
            # indices2 = indices1 + 1 + np.random.choice(data_length - 1,
            #                                            size=n_cannot_link)
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
