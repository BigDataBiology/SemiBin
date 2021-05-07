import torch
from torch.nn import Linear, ReLU, LeakyReLU
from torch import nn
from torch.utils.data import Dataset, DataLoader
import os
from .utils import cal_num_bins
import pandas as pd
from torch.optim import lr_scheduler
from tqdm import tqdm


class Semi_encoding_multiple(torch.nn.Module):
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


class Semi_encoding_single(torch.nn.Module):
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

    def __init__(self, embedding1_lists, embedding2_lists, label_lists):
        self.embedding1_lists = embedding1_lists
        self.embedding2_lists = embedding2_lists
        self.label_lists = label_lists

    def __getitem__(self, item):
        return self.embedding1_lists[item], self.embedding2_lists[item], self.label_lists[item]

    def __len__(self):
        return len(self.embedding1_lists)


class unsupervised_feature_Dataset(Dataset):
    def __init__(self, embedding1_lists, embedding2_lists):
        self.embedding1_lists = embedding1_lists
        self.embedding2_lists = embedding2_lists

    def __getitem__(self, item):
        return self.embedding1_lists[item], self.embedding2_lists[item]

    def __len__(self):
        return len(self.embedding1_lists)


def train(out, contig_fasta, binned_short, logger, data, data_split, cannot_link, is_combined=True,
          batchsize=2048, epoches=20, device=None,num_process = 8):
    contig_output = os.path.join(out, os.path.split(contig_fasta)[1] + '.frag')
    hmm_output = os.path.join(out, os.path.split(contig_fasta)[1] + '.hmmout')
    seed_output = os.path.join(out, os.path.split(contig_fasta)[1] + '.seed')


    cal_num_bins(
        contig_fasta,
        contig_output,
        hmm_output,
        seed_output,
        binned_short,
        num_process)



    logger.info('Generate training data:')

    # #generate two inputs
    train_input_1 = []
    train_input_2 = []
    train_labels = []

    # cannot link
    cannot_link = pd.read_csv(cannot_link, sep=',', header=None).values

    namelist = data.index.tolist()
    mapObj = dict(zip(namelist, range(len(namelist))))
    train_data = data.values
    train_data_must_link = data_split.values

    if not is_combined:
        train_data_input = train_data[:, 0:136]
        train_data_split_input = train_data_must_link
    else:
        train_data_input = train_data
        train_data_split_input = train_data_must_link

    # can not link from contig annotation
    for link in cannot_link:
        train_input_1.append(train_data_input[mapObj[str(link[0])]])
        train_input_2.append(train_data_input[mapObj[str(link[1])]])
        train_labels.append(0)

    # cannot link from bin seed
    if os.path.exists(seed_output):
        seed = open(seed_output).read().split('\n')
        seed = [contig for contig in seed if contig != '']
        for i in range(len(seed)):
            for j in range(i + 1, len(seed)):
                train_input_1.append(train_data_input[mapObj[str(seed[i])]])
                train_input_2.append(train_data_input[mapObj[str(seed[j])]])
                train_labels.append(0)

    # must link from breaking up
    for i in range(0, len(train_data_must_link), 2):
        train_input_1.append(train_data_split_input[i])
        train_input_2.append(train_data_split_input[i + 1])
        train_labels.append(1)

    logger.info('Number of must link pair:{}'.format(train_labels.count(1)))
    logger.info('Number of can not link pair:{}'.format(train_labels.count(0)))

    logger.info('Training model...')

    dataset = feature_Dataset(train_input_1, train_input_2, train_labels)
    train_loader = DataLoader(
        dataset=dataset,
        batch_size=batchsize,
        shuffle=True,
        num_workers=0)

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
        num_workers=16)

    if not is_combined:
        model = Semi_encoding_single(train_data_input.shape[1]).to(device)
    else:
        model = Semi_encoding_multiple(train_data_input.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)

    for epoch in tqdm(range(epoches)):
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
    torch.save(model, os.path.join(out, 'model.h5'))

    return model