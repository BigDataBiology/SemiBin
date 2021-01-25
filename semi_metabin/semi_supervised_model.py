import torch
from torch.nn import Linear,ReLU,LeakyReLU
from torch import nn
from torch.utils.data import Dataset,DataLoader

class Semi_encoding_multiple(torch.nn.Module):
    def __init__(self,num):
        super(Semi_encoding_multiple,self).__init__()
        self.encoder1 = torch.nn.Sequential(
            Linear(num,512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512,100),
        )


        self.decoder1 = torch.nn.Sequential(
            Linear(100,512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512,num),
            nn.Sigmoid(),
        )

    def forward(self,input1,input2):
        return self.encoder1(input1) , self.encoder1(input2)

    def decoder(self,input1,input2):
        return self.decoder1(input1) , self.decoder1(input2)


    def embedding(self,input):
        return self.encoder1(input)

class Semi_encoding_single(torch.nn.Module):
    def __init__(self,num):
        super(Semi_encoding_single,self).__init__()
        self.encoder1 = torch.nn.Sequential(
            Linear(num,512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512,100),
            #nn.Sigmoid(),
        )


        self.decoder1 = torch.nn.Sequential(
            Linear(100,512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512, 512),
            nn.BatchNorm1d(512),
            LeakyReLU(),
            nn.Dropout(0.2),
            Linear(512,num),
            nn.Softmax(dim=1),
        )



    def forward(self,input1,input2):
        return self.encoder1(input1) , self.encoder1(input2)

    def decoder(self,input1,input2):
        return self.decoder1(input1) , self.decoder1(input2)


    def embedding(self,input):
        return self.encoder1(input)



def loss_function(embedding1,embedding2,label,raw_x_1,raw_x_2,decoder_x_1,decoder_x_2,is_label = True):
    relu = torch.nn.ReLU()
    mse_loss = torch.nn.MSELoss()
    d = torch.norm(embedding1-embedding2,p=2,dim=1)
    sqaure_pred = torch.square(d)
    margin_square = torch.square(relu(1-d))

    if is_label:
        supervised_loss = torch.mean(label * sqaure_pred + (1 - label) * margin_square)
        unsupervised_loss = 0.5 * mse_loss(decoder_x_1,raw_x_1) + 0.5 * mse_loss(decoder_x_2,raw_x_2)
        loss = supervised_loss + unsupervised_loss
        return loss,supervised_loss,unsupervised_loss

    else:
        unsupervised_loss = 0.5 * mse_loss(decoder_x_1, raw_x_1) + 0.5 * mse_loss(decoder_x_2, raw_x_2)
        return unsupervised_loss

class feature_Dataset(Dataset):

    def __init__(self, embedding1_lists, embedding2_lists,label_lists):
        self.embedding1_lists = embedding1_lists
        self.embedding2_lists = embedding2_lists
        self.label_lists = label_lists

    def __getitem__(self, item):
        return self.embedding1_lists[item], self.embedding2_lists[item],self.label_lists[item]

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