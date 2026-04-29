import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import torch
import torch.nn as nn
import torch.optim as optim

from src.data_gen import DataGen
from src.train import TrainData
from src.models.PyomoLayers import PyomoLayers

from src.utils_c_hzn import (
    policy,
    load_data,
    format_data,
    build_train_val_dataset,
    build_codesign_dataset,
    train_dac_surrogate,
    formulate_OMLT,
    formulate_model,
    # run_codesign,
    # vis_codesign,
    test_model,
    vis_test
)

PATH = '/Users/dwald/Documents/DOL/DOL-DAC/'

###----------------------------- LOADING THE DAC TRAINING DATASET -----------------------------###
filename = 'XY_13D_N11_30k_K_4_N_120k_updated.csv'
X, y = load_data(PATH, filename)

###----------------------------- FORMATTING THE DAC TRAINING DATASET -----------------------------###
X_train_norm_torch, Y_train_torch, X_test_norm_torch, Y_test_torch, scaler = format_data(X,y)

###----------------------------- DAC TRAINING DATASET IN TORCH FORMAT -----------------------------###
batch_size = 100
dataloader_train, dataloader_val = build_train_val_dataset(X_train_norm_torch, Y_train_torch, X_test_norm_torch, Y_test_torch, batch_size=batch_size)

###----------------------------- TRAIN THE DAC SURROGATE MODEL -----------------------------###
in_dim = X.shape[1]
hidden_dim = in_dim*3
out_dim = y.shape[1]
dac = train_dac_surrogate(dataloader_train, dataloader_val, in_dim, hidden_dim, out_dim, PATH, save=False)

###----------------------------- LOAD THE DAC SURROGATE MODEL INTO OMLT -----------------------------###
Nh = 2
lb = np.min(X_train_norm_torch.detach().numpy(), axis=0).astype(np.float64)#; lb[0] = -0.2; lb[1] = -0.2
ub = np.max(X_train_norm_torch.detach().numpy(), axis=0).astype(np.float64)#; ub[0] = 1.2; ub[1] = 1.2
input_bounds = list(zip(lb, ub))
# formulation = formulate_OMLT(dac, input_bounds)
formulations = [formulate_OMLT(dac, input_bounds) for _ in range(Nh)]

###----------------------------- CREATE THE PYOMOLAYER -----------------------------###
# obj='track_demand'
obj='minimize_cost'
# obj='minimize_net'
# obj='combined'
# obj='min_demand'
# obj='max_productivity'
m = formulate_model(formulations, Nh, obj=obj)

###----------------------------- RUN THE CONTROL CO-DESIGN PROCESS -----------------------------###
batch_size_codesign = 1#876
filename_codesign = "data/weather_data_WY.csv"
dataloader_codesign, ref_all = build_codesign_dataset(filename_codesign, batch_size_codesign, scaler, Nh, obj=obj)

# epochs = 1#10
# # lr = 1e-3 # for 'track demand'
# # lr = 1e-4 # for 'minimize cost'
# lr = 0.0 #8e-5 # for 'combined'
# # lr = 8e-5 # for 'min demand'
# # lr = 5e-5 # for 'max poductivity

# if (obj == 'track_demand') or (obj == 'minimize_cost'):
#     in_dim = 3
# else:
#     in_dim = 2
# policy_net = policy(in_dim=in_dim, out_dim=int(len(Pyomolayer.concrete_model.ctrl_vars)/Nh), hidden_dim=100).float()
# history, history_u, policy_net = run_codesign(Pyomolayer, policy_net, dataloader_codesign, scaler, param_vec_initial, Nh, obj=obj,
#                                                 batch_size=batch_size_codesign, lr=lr, epochs=epochs,
#                                                 batch_size_u=batch_size_u, lr_u=lr_u, epochs_u=epochs_u)

###----------------------------- VISUALIZE THE TRAINING RESULTS -----------------------------###
# vis_codesign(PATH, history, history_u, save_csv=True)

###----------------------------- TEST THE SYSTEM WITH REAL ENVIRONMENTAL CONDITIONS -----------------------------###
filename_test = "data/weather_data_WY.csv"
trajectories = test_model(dac, m, filename_test, scaler, Nh, obj, ref_all)

###----------------------------- VISUALIZE THE TEST RESULTS -----------------------------###
vis_test(PATH, trajectories, show=True, save=True)