

'''
Installation:
module load python/3.8-anaconda3
conda create --name light3 python=3.8
conda activate light3
pip install pytorch-lightning
pip install torchvision
pip install scikit-learn
'''
import torch
from torch import nn
from torch.nn import functional as F
from torch.utils.data import DataLoader
from torch.utils.data import random_split
from torchvision import transforms
import pytorch_lightning as pl
from pytorch_lightning.plugins import DDPPlugin

# imports from stemdl
import time,sys, os, math, glob,argparse, yaml
import torch.backends.cudnn as cudnn
import torch.multiprocessing as mp
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data.distributed
from torch.utils.tensorboard import SummaryWriter
from torchvision import datasets, transforms, models
from tqdm import tqdm
from sklearn.metrics import f1_score
import torch.nn as nn
import numpy as np

from pathlib import Path
from torch.utils.data import Dataset
from torchvision import datasets
from torchvision.transforms import ToTensor

# Custom dataset class
class NPZDataset(Dataset):
    def __init__(self, npz_root):
        self.files = glob.glob(npz_root + "/*.npz")
        #self.files = files1[0:63]

    def __getitem__(self, index):
        sample = np.load(self.files[index])
        x = torch.from_numpy(sample["data"])
        y = sample["label"][0]
        return (x, y)

    def __len__(self):
        return len(self.files)

# StemdlModel
class StemdlModel(pl.LightningModule):
    def __init__(self):
        super().__init__()
        self.input_size = 128
        self.num_classes = 231
        self.model_name = "resnet18"
        self.model = models.resnet18(pretrained=False) 
        self.num_ftrs = self.model.fc.in_features
        self.model.fc = nn.Linear(self.num_ftrs, self.num_classes)
        self.params_to_update = self.model.parameters()
        self.feature_extract = False

    # forward step
    def forward(self, x):
        embedding = self.model(x)
        return embedding

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        return optimizer

    def training_step(self, train_batch, batch_idx):
        x, y = train_batch
        x_hat = self.model(x)
        y = F.one_hot(y, num_classes=231).float()
        loss = F.mse_loss(x_hat, y)
        self.log('train_loss', loss)
        return loss

    def validation_step(self, val_batch, batch_idx):
        x, y = val_batch
        x_hat = self.model(x)
        y = F.one_hot(y, num_classes=231).float()
        loss = F.mse_loss(x_hat, y)
        self.log('train_loss', loss)
        return loss
    
    def test_step(self, test_batch, batch_idx):
        x, y = test_batch
        x_hat = self.model(x)
        y = F.one_hot(y, num_classes=231).float()
        loss = F.mse_loss(x_hat, y)
        self.log('test_loss', loss)
        return loss

    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x, y = batch
        y_hat = self.model(x)
        return y_hat

# Running the code: 
# python stemdl_light.py 
#
def main():

    # Read YAML file
    with open("stemdlConfig.yaml", 'r') as stream:
        config = yaml.safe_load(stream)

    # Datasets: training (138717 files), validation (20000 files),
    # testing (20000 files), prediction (8438 files), 197kbytes each
    train_dataset = NPZDataset(os.path.expanduser(config['train_dir']))
    val_dataset = NPZDataset(os.path.expanduser(config['val_dir']))
    test_dataset = NPZDataset(os.path.expanduser(config['test_dir']))
    predict_dataset = NPZDataset(os.path.expanduser(config['inference_dir']))
    
    # Data loaders
    train_loader = DataLoader(dataset=train_dataset,batch_size=int(config['batchsize']),num_workers=4)

    val_loader = DataLoader(dataset=val_dataset,batch_size=int(config['batchsize']),num_workers=4)

    test_loader = DataLoader(dataset=test_dataset,batch_size=int(config['batchsize']),num_workers=4)

    predict_loader = DataLoader(dataset=predict_dataset,batch_size=int(config['batchsize']),num_workers=4)

    # model
    model = StemdlModel()

    # training
    trainer = pl.Trainer(gpus=int(config['gpu']), num_nodes=int(config['nodes']), precision=16, strategy="ddp", max_epochs=int(config['epochs']))
    start = time.time()
    trainer.fit(model, train_loader, val_loader)
    diff = time.time() - start
    elapsedTime =  f"{diff:.2f}"

    log_file = os.path.expanduser(config['log_file'])
    with open(log_file, "a") as logfile:
        logfile.write(f"Stemdl, resnet={config['resnet']}, epochs={config['epochs']}, bs={config['batchsize']}, nodes={config['nodes']}, gpu={config['gpu']}, time={elapsedTime}\n")
    
    #testing
    trainer.test(model, test_loader)

    #inference
    predictions = trainer.predict(model, dataloaders=predict_loader)
    
if __name__ == "__main__":
        main()
