import os
import pickle
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from timm.utils import ModelEmaV3
from torch.utils.data import DataLoader
from tqdm import tqdm

from Diffusion import set_seed, DDPM_Scheduler, UNET


def diffusion(batch_size: int = 16,
              num_time_steps: int = 1000,
              num_epochs: int = 150000,
              save_epoch: int = 10000,
              add_method= None,
              add_epochs: int = 0,
              seed: int = 42,
              ema_decay: float = 0.9999,
              lr=1e-5,
              train_type='case',  # case or ctrl
              device=None,
              formatted_time='0',
              checkpoint_path: str = None):
    set_seed(seed)
    data_o_case, data_o_ctrl, taxa_list = load_sample_pickle_data("checkpoints/raw_data.pkl")
    if train_type == "case":
        data_o = data_o_case.astype(float)
    else:
        data_o = data_o_ctrl.astype(float)
    data_reshape_matrix = ImageCoding(data_o, 32, 24)
    train_loader = DataLoader(data_reshape_matrix, batch_size=batch_size, shuffle=True, drop_last=True, num_workers=0)
    scheduler = DDPM_Scheduler(num_time_steps=num_time_steps)
    model = UNET().cuda()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    ema = ModelEmaV3(model, decay=ema_decay)
    if checkpoint_path is not None:
        checkpoint = torch.load(checkpoint_path)
        model.load_state_dict(checkpoint['weights'])
        ema.load_state_dict(checkpoint['ema'])
        optimizer.load_state_dict(checkpoint['optimizer'])
    criterion = nn.MSELoss(reduction='mean')

    for epoch in range(add_epochs, add_epochs + num_epochs):
        with tqdm(train_loader, dynamic_ncols=False) as tqdmDataLoader:
            for otu in tqdmDataLoader:
                x = otu.float().to(device)
                t = torch.randint(0, num_time_steps, (batch_size,))
                e = torch.randn_like(x, requires_grad=False).to(device)
                a = scheduler.alpha[t].view(batch_size, 1, 1, 1).to(device)
                x = (torch.sqrt(a) * x) + (torch.sqrt(1 - a) * e)
                output = model(x, t)
                optimizer.zero_grad()
                loss = criterion(output, e)
                loss.backward()
                optimizer.step()
                ema.update(model)
                tqdmDataLoader.set_postfix(ordered_dict={
                    "epoch": str(epoch + 1) + "/" + str(add_epochs + num_epochs),
                    "loss: ": loss.item(),
                    "x.shape: ": x.shape,
                    "LR": optimizer.state_dict()['param_groups'][0]["lr"],
                })
        if (epoch + 1) % save_epoch == 0:
            checkpoint = {
                'weights': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'ema': ema.state_dict()
            }
            torch.save(checkpoint,
                       os.path.join('.//checkpoints/epoch_' + str(epoch + 1) + "_step" + str(num_time_steps) +
                                    "_" + str(train_type) + '_' + str(add_method) + "_" + str(formatted_time) + ".pt"))


def reverse_diffusion(checkpoint_path: str = None,
                     num_time_steps: int = 1000,
                     ema_decay: float = 0.9999,
                     train_type='case',
                     generate_num=1000,
                     formatted_time='0',
                     num_epochs=0,
                     add_method=None,
                     device='cuda:0'):
    checkpoint = torch.load(checkpoint_path, weights_only=True)
    model = UNET().cuda()
    model.load_state_dict(checkpoint['weights'])
    ema = ModelEmaV3(model, decay=ema_decay)
    ema.load_state_dict(checkpoint['ema'])
    scheduler = DDPM_Scheduler(num_time_steps=num_time_steps)
    results = []

    with torch.no_grad():
        model = ema.module.eval()
        for i in range(generate_num):
            z = torch.randn(1, 1, 32, 24)
            for t in reversed(range(1, num_time_steps)):
                t = [t]
                temp = (scheduler.beta[t] / (
                        (torch.sqrt(1 - scheduler.alpha[t])) * (torch.sqrt(1 - scheduler.beta[t]))))
                z = (1 / (torch.sqrt(1 - scheduler.beta[t]))) * z - (temp * model(z.cuda(), t).cpu())
                e = torch.randn(1, 1, 32, 24)
                z = z + (e * torch.sqrt(scheduler.beta[t]))
            temp = scheduler.beta[0] / ((torch.sqrt(1 - scheduler.alpha[0])) * (torch.sqrt(1 - scheduler.beta[0])))
            x = (1 / (torch.sqrt(1 - scheduler.beta[0]))) * z - (temp * model(z.cuda(), [0]).cpu())
            x_flattened = x.view(1, -1)
            results.append(x_flattened)
            print(i + 1)
    title_filename = str(".//checkpoints//" + train_type + "_title.csv")
    title = pd.read_csv(title_filename, header=None).iloc[0]
    results = torch.cat(results, dim=0)[:, :title.shape[0]].numpy()
    df = pd.DataFrame(results, columns=title)
    df.to_csv('.//data//epoch_' + str(num_epochs) + '_' + str(train_type) + '_' + str(add_method) + '_' + str(
        formatted_time) + '.csv', index=False)

def load_sample_pickle_data(filename=".//checkpoints//raw_data.pkl"):
    # Load raw dataset
    raw_data = pickle.load(open(filename, 'rb'))
    dataset = raw_data.iloc[:, 1:].values
    labels = raw_data["group"].values
    taxa_list = raw_data.columns[1:]

    data_o_case = dataset[labels == 'case']
    data_o_ctrl = dataset[labels == 'ctrl']

    return data_o_case, data_o_ctrl, taxa_list

def ImageCoding(matrix, m, n):
    target_dim = m * n
    if matrix.shape[1] > target_dim:
        raise ValueError("The number of columns in the input matrix cannot exceed the target dimension:" + str(m) + "*" + str(n))
    padded_matrix = np.zeros((matrix.shape[0], target_dim))
    padded_matrix[:, :matrix.shape[1]] = matrix  # 填充原始数据
    tensor = padded_matrix.reshape(matrix.shape[0], 1, m, n)
    return tensor