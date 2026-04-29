import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as tm
import random

import warnings

from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

import torch
import torch.nn as nn

import pyomo.environ as pyo
from pyomolayers import PyomoOptLayer

import torch.onnx
from omlt.io import load_onnx_neural_network_with_bounds, write_onnx_model_with_bounds
import tempfile
from omlt import OmltBlock
from omlt.neuralnet import FullSpaceNNFormulation


class sys(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim):
        super().__init__()

        self.in_layer = nn.Linear(in_dim, hidden_dim)
        self.hidden_layer = nn.Linear(hidden_dim, hidden_dim)
        self.out_layer = nn.Linear(hidden_dim, out_dim)

        self.act = nn.Sigmoid()

    def forward(self, x):

        h = self.act(self.in_layer(x))
        h = self.act(self.hidden_layer(h))
        y = self.out_layer(h)

        return y
    
class policy(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim):
        super().__init__()

        self.in_layer = nn.Linear(in_dim, hidden_dim)
        self.hidden_layer = nn.Linear(hidden_dim, hidden_dim)
        self.out_layer = nn.Linear(hidden_dim, out_dim)

        self.act = nn.Sigmoid()

    def forward(self, x):

        h = self.act(self.in_layer(x))
        h = self.act(self.hidden_layer(h))
        y = self.out_layer(h)

        return y
    
def build_reference(obj):

    if obj == 'minimize_net':
        ### hourly carbon intensity data
        # dCO2_df = pd.read_csv('/Users/dwald/Documents/DOL/DOL-DAC/data/Cambium24_MidCase_hourly_p61_2025.csv') # for TX
        dCO2_df = pd.read_csv('/Users/dwald/Documents/DOL/DOL-DAC/data/Cambium24_MidCase_hourly_p22_2025.csv') # for WY
        intensity_all = dCO2_df[dCO2_df.columns[50]].iloc[2:].to_numpy().astype(float).reshape(-1,1) # kgCO2/MWh
        intensity_all = intensity_all * (1./1000.) # tCO2/MWh

        return intensity_all

    if obj == 'minimize_cost':
        ### hourly carbon intensity data
        # dCO2_df = pd.read_csv('/Users/dwald/Documents/DOL/DOL-DAC/data/Cambium24_MidCase_hourly_p61_2025.csv') # for TX
        dCO2_df = pd.read_csv('/Users/dwald/Documents/DOL/DOL-DAC/data/Cambium24_MidCase_hourly_p22_2025.csv') # for WY
        intensity_all = dCO2_df[dCO2_df.columns[50]].iloc[2:].to_numpy().astype(float).reshape(-1,1) # kgCO2/MWh
        intensity_all = intensity_all * (1./1000.) # tCO2/MWh

        ### hourly cost of energy prices
        # cost_df = pd.read_excel('/Users/dwald/Documents/DOL/DOL-DAC/data/TX_WY_2025.xlsx', sheet_name='Atoka TX 2025') # for TX
        cost_df = pd.read_excel('/Users/dwald/Documents/DOL/DOL-DAC/data/TX_WY_2025.xlsx', sheet_name='Madison WY 2025') # for WY
        cost_all = cost_df['$/MWh'].to_numpy().astype(float).reshape(-1,1) # $/MWh
        # cost_all = cost_all * (1/100.)

        intensity_cost_all = np.hstack((intensity_all, cost_all))

        return intensity_cost_all
        
    if obj == 'track_demand':
        val = 0.25
        ref_all = np.zeros((8760,1))
        for i in range(len(ref_all)):
            if i%730 == 0.0:
                val = 0.25 + random.uniform(-0.05, 0.05)
            ref_all[i,:] = val
        
        return ref_all

def load_data(PATH, filename):

    df = pd.read_csv(os.path.join(PATH,'data',filename))

    # convert total cycle time from seconds to hours
    df['t_cycle_h'] = df['t_cycle_s'] / 60. / 60. # convert to hours

    # df['RH_in_frac'] = df['RH_in_frac'] / 100.
    
    X_cols = ['T_ambC_C', 'RH_in_frac','qCO2_start_molkg', 'qH2O_start_molkg', 'u_feed_m_s', 'p_vac_bar', 'T_reg_C','t_heat_s', 't_des_s', 't_ads_s', 'width_m', 'height_m', 'N_cell']
    # y_cols = ['qCO2_postAds_molkg','qCO2_postDes_molkg', 'qH2O_postAds_molkg', 'qH2O_postDes_molkg','E_cycle_MWh_cycle_collector']
    y_cols = ['qCO2_postDes_molkg', 'qH2O_postDes_molkg','m_CO2 (t-CO2/cycle)','E_cycle_MWh_cycle','t_cycle_h']
    X, y = df[X_cols].to_numpy(), df[y_cols].to_numpy()

    return X, y

def format_data(X, y):

    X_train, X_test, Y_train, Y_test = train_test_split(X,y,train_size=0.95)
    print(f"Surrogate model training data shape: X_train: {X_train.shape}, Y_train: {Y_train.shape}, X_test: {X_test.shape}, Y_test: {Y_test.shape}")
    
    # Normalize input data
    scaler = MinMaxScaler()
    X_train_norm = scaler.fit_transform(X_train)
    X_test_norm = scaler.transform(X_test)

    # Normalize the state variables only in the output dataset (for feedback)
    Y_train[:,0:2] = scaler.transform(np.hstack((np.zeros((Y_train.shape[0],2)),Y_train[:,0:2],np.zeros((Y_train.shape[0],6)),np.zeros((Y_train.shape[0],3)))))[:,2:4]
    Y_test[:,0:2] = scaler.transform(np.hstack((np.zeros((Y_test.shape[0],2)),Y_test[:,0:2],np.zeros((Y_test.shape[0],6)),np.zeros((Y_test.shape[0],3)))))[:,2:4]

    # Build the torch dataset
    X_train_norm_torch = torch.tensor(X_train_norm, dtype=torch.float32)
    Y_train_torch = torch.tensor(Y_train, dtype=torch.float32)
    X_test_norm_torch = torch.tensor(X_test_norm, dtype=torch.float32)
    Y_test_torch = torch.tensor(Y_test, dtype=torch.float32)
    # X_train_norm_torch.shape, Y_train_torch.shape, X_test_norm_torch.shape, Y_test_torch.shape
    return X_train_norm_torch, Y_train_torch, X_test_norm_torch, Y_test_torch, scaler

def build_train_val_dataset(X_train, Y_train, X_test, Y_test, batch_size):

    dataset_train = torch.utils.data.TensorDataset(X_train, Y_train)
    dataloader_train = torch.utils.data.DataLoader(dataset_train, batch_size=batch_size, shuffle=False)

    dataset_val = torch.utils.data.TensorDataset(X_test, Y_test)
    dataloader_val = torch.utils.data.DataLoader(dataset_val, batch_size=1, shuffle=False)

    return dataloader_train, dataloader_val

def train_dac_surrogate(dataloader_train, dataloader_val, in_dim, hidden_dim, out_dim, PATH, save=False):
    if save is True:
        dac = sys(in_dim, hidden_dim, out_dim).float()

        learn_rate = 1e-4 #1e-3
        epochs = 5000 #1500

        opt = torch.optim.Adam(dac.parameters(), lr=learn_rate)

        t1 = tm.time()
        loss_epoch = []
        val_loss_epoch = []
        for epoch in range(epochs):
            
            dac.train()
            loss_batch_train = []
            for batch_idx, (x, y) in enumerate(dataloader_train):
                opt.zero_grad()
                y_hat = dac(x)
                loss = (y_hat - y).pow(2).mean()
                loss.backward()
                opt.step()
                loss_batch_train.append(loss.item())

            dac.eval()
            loss_batch_val = []
            with torch.no_grad():
                for batch_idx, (x_val, y_val) in enumerate(dataloader_val):
                    y_hat_val = dac(x_val)
                    val_loss = (y_hat_val - y_val).pow(2).mean()
                    loss_batch_val.append(val_loss.item())

            if epoch%100 == 0:
                print(f"Epoch {epoch} training loss: {np.mean(loss_batch_train):.3e}, validation loss: {np.mean(loss_batch_val):.3e}")
            loss_epoch.append(np.mean(loss_batch_train))
            val_loss_epoch.append(np.mean(loss_batch_val))

        t2 = tm.time()
        print(f"Total training time: {t2-t1:.3f} sec")
        print(f"Final train data mean loss: {loss_epoch[-1]}")
        print(f"Final validation data mean loss: {val_loss_epoch[-1]}")

        plt.figure()
        plt.plot(loss_epoch, label="Training")
        plt.plot(val_loss_epoch, label="Validation")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.yscale('log')
        plt.legend()
        plt.show()
        
        torch.save(
            dac.state_dict(),
            os.path.join(PATH,'src','models','dac_surrogate_new')
            )

        return dac.eval()

    if save is False:
        dac = sys(in_dim, hidden_dim, out_dim)
        dac.load_state_dict(torch.load(os.path.join(PATH,'src','models','dac_surrogate_new'), weights_only=True))
        return dac.eval().float()
    

def formulate_OMLT(dac, input_bounds):

    # model input used for exporting
    x = torch.randn(1, len(input_bounds), requires_grad=True, dtype=torch.float32)
    pytorch_model = None
    with tempfile.NamedTemporaryFile(suffix=".onnx", delete=False) as f:
        torch.onnx.export(
            dac,
            x,
            f.name,
            input_names=["input"],
            output_names=["output"],
            dynamic_axes={"input": {0: "batch_size"}, "output": {0: "batch_size"}},
        )
        write_onnx_model_with_bounds(f.name, None, input_bounds)
        print(f"Wrote PyTorch model to {f.name}")
        pytorch_model = f.name

    network_definition = load_onnx_neural_network_with_bounds(pytorch_model)

    formulation = FullSpaceNNFormulation(network_definition)
    
    return formulation

def build_codesign_dataset(filename_cd, batch_size, scaler, Nh, obj):

    # Nh = 12
    df = pd.read_csv(filename_cd)
    df["Relative Humidity (%)"] = df["Relative Humidity (%)"] / 100 # convert from percentage to decimal

    lower_bounds = {"Temperature (degC)" : 5.0, "Relative Humidity (%)" : 0.1}
    upper_bounds = {"Temperature (degC)" : 40.0, "Relative Humidity (%)" : 0.9}
    lower = pd.Series(lower_bounds)
    upper = pd.Series(upper_bounds)
    df = df.clip(lower=lower, upper=upper, axis=1)

    X = df[["Temperature (degC)","Relative Humidity (%)"]].to_numpy()
    y_dummy = np.zeros_like(X)

    # X_codesign, _, Y_codesign, _ = train_test_split(X,y_dummy, train_size=0.9,shuffle=True,random_state=101)
    X_codesign, Y_codesign = X, y_dummy

    X_codesign_norm_all = scaler.transform(np.hstack((X_codesign, np.zeros((X_codesign.shape[0], 11)))))
    X_codesign_norm = X_codesign_norm_all[:,:2]
    
    if obj == 'minimize_net':
        intensity_all = build_reference(obj)
        X_codesign_norm = np.hstack((X_codesign_norm, intensity_all))
    if obj == 'minimize_cost':
        intensity_cost_all = build_reference(obj)
        X_codesign_norm = np.hstack((X_codesign_norm, intensity_cost_all))
    if obj == 'track_demand':
        ref_all = build_reference(obj)
        X_codesign_norm = np.hstack((X_codesign_norm, ref_all))
    else:
        ref_all = None
    
    X_codesign_norm_hzn = np.zeros((X_codesign_norm.shape[0], Nh, X_codesign_norm.shape[1]))
    for i in range(X.shape[0]-Nh):
        X_codesign_norm_hzn[i, :, :] = X_codesign_norm[i:i+Nh, :]
    Y_codesign_hzn = np.zeros_like(X_codesign_norm_hzn)
    
    print(f"Co-design training data shape: X: {X_codesign_norm_hzn.shape} (no y-data needed)")

    X_codesign_norm_torch = torch.tensor(X_codesign_norm_hzn, dtype=torch.float32)
    Y_codesign_torch = torch.tensor(Y_codesign_hzn, dtype=torch.float32)

    dataset_codesign = torch.utils.data.TensorDataset(X_codesign_norm_torch, Y_codesign_torch)
    dataloader_codesign = torch.utils.data.DataLoader(dataset_codesign, batch_size=batch_size, shuffle=False)

    return dataloader_codesign, ref_all

def initial_time(model):
    return model.Ts[0] == 0.0
def time_transition(model):
    return model.Ts[1] == model.nn1.outputs[4] + model.Ts[0]

def initial_env(model, i):
    return model.env_vars[0, i] == model.env_vars0[i]
def calc_env_cond(model, i):
    return model.env_vars[1, i] == model.beta1[i] * model.Ts[1]**2 + model.beta2[i] * model.Ts[1] + model.beta3[i]
def calc_reference(model, k):
    return model.ref[k] == ((model.K-model.C) / (1 + pyo.exp(model.k*(model.Ts[k] - model.Ts_0)))) + model.C

def initial_state(model, k):
    return model.state_vars[0,k] == model.state_0_vars[k]
def state_transition(model, k):
    return model.state_vars[1,k] == model.nn1.outputs[k]

def output1_lower_bound(model, k):
    return model.nn1.outputs[k] >= 0.0
def output2_lower_bound(model, k):
    return model.nn2.outputs[k] >= 0.0

def input1_rule_lower(model, k):
    return model.nn1.inputs[k] >= 0.0
def input2_rule_lower(model, k):
    return model.nn2.inputs[k] >= 0.0

def input1_rule_upper(model, k):
    return model.nn1.inputs[k] <= 1.0
def input2_rule_upper(model, k):
    return model.nn2.inputs[k] <= 1.0

def formulate_model(formulations, Nh, obj):

    # create pyomo model
    m = pyo.ConcreteModel()

    # create an OMLT block for the neural network and build its formulation
    m.nn1 = OmltBlock()
    m.nn1.build_formulation(formulations[0])
    m.nn2 = OmltBlock()
    m.nn2.build_formulation(formulations[1])

    m.Ts = pyo.Var(range(Nh))
    m.beta1 = pyo.Var(range(2))
    m.beta2 = pyo.Var(range(2))
    m.beta3 = pyo.Var(range(2))
    m.env_vars0 = pyo.Var(range(2))
    
    # separate nn inputs int environment variables, control variables, and parameters
    N_envs = 2
    N_states = 2
    N_ctrls = 6
    N_params = 3
    m.env_vars   = pyo.Var(range(Nh), range(N_envs))
    m.state_vars   = pyo.Var(range(Nh), range(N_states))
    m.ctrl_vars  = pyo.Var(range(Nh), range(N_ctrls))
    m.param_vars = pyo.Var(range(Nh), range(N_params))
    if obj == 'track_demand':
        m.ref = pyo.Var(range(Nh))
        m.K = pyo.Var()
        m.C = pyo.Var()
        m.Ts_0 = pyo.Var()
        m.k = pyo.Var()
    if obj == 'minimize_net':
        m.intensity = pyo.Var(range(Nh))
    if obj == 'minimize_cost':
        m.intensity = pyo.Var(range(Nh))
        m.cost = pyo.Var(range(Nh))
    x1 = m.nn1.inputs
    x2 = m.nn2.inputs
    # envs
    env_ids = [0,1]
    for i,id in enumerate(env_ids):
        m.add_component(f"map_env1_{i}", pyo.Constraint(expr=x1[id] == m.env_vars[0,i]))
        m.add_component(f"map_env2_{i}", pyo.Constraint(expr=x2[id] == m.env_vars[1,i]))
    # states
    state_ids = [2,3]
    for i,id in enumerate(state_ids):
        m.add_component(f"map_state1_{i}", pyo.Constraint(expr=x1[id] == m.state_vars[0,i]))
        m.add_component(f"map_state2_{i}", pyo.Constraint(expr=x2[id] == m.state_vars[1,i]))
    # ctrls
    ctrl_ids = [4,5,6,7,8,9]
    for j,id in enumerate(ctrl_ids):
        m.add_component(f"map_ctrl1_{j}", pyo.Constraint(expr=x1[id] == m.ctrl_vars[0,j]))
        m.add_component(f"map_ctrl2_{j}", pyo.Constraint(expr=x2[id] == m.ctrl_vars[1,j]))
    # params
    param_ids = [10,11,12]
    for k, id in enumerate(param_ids):
        m.add_component(f"map_param1_{k}", pyo.Constraint(expr=x1[id] == m.param_vars[0,k]))
        m.add_component(f"map_param2_{k}", pyo.Constraint(expr=x2[id] == m.param_vars[1,k]))
    m.state_0_vars = pyo.Var(range(N_states))

    # iniitalize the environment variables and parameters (parameters in Pyomolayer model)
    # beta1_init = np.random.rand(2)
    # for j in range(len(beta1_init)):
    #     m.beta1[j].fix(beta1_init[j])
    # beta2_init = np.random.rand(2)
    # for j in range(len(beta2_init)):
    #     m.beta2[j].fix(beta2_init[j])
    # beta3_init = np.random.rand(2)
    # for j in range(len(beta3_init)):
    #     m.beta3[j].fix(beta3_init[j])
    # env_var0_init = np.random.rand(2)
    # for j in range(len(env_var0_init)):
    #     m.env_vars0[j].fix(env_var0_init[j])
    # env_var_init = np.random.rand(Nh, len(env_ids))
    # for i in range(Nh):
    #     for j in range(len(env_ids)):
    #         m.env_vars[i,j].fix(env_var_init[i,j])
    # state_0_var_init = np.random.rand(len(state_ids))
    # for j in range(len(state_ids)):
    #     m.state_0_vars[j].fix(state_0_var_init[j])
    param_var_init = np.array([[0.5,0.5,0.5],
                               [0.5,0.5,0.5]])
    for i in range(Nh):
        for j in range(len(param_ids)):
            m.param_vars[i,j].fix(param_var_init[i,j])
    # if obj == 'track_demand':
    #     K_init = np.random.rand(1)
    #     C_init = np.random.rand(1)
    #     Ts_0_init = np.random.rand(1)
    #     k_init = np.random.rand(1)
    #     m.K.fix(K_init)
    #     m.C.fix(C_init)
    #     m.Ts_0.fix(Ts_0_init)
    #     m.k.fix(k_init)
    # if obj == 'minimize_cost':
    #     weight_init = np.random.rand(2)
    #     for j in range(len(weight_init)):
    #             m.weight[j].fix(weight_init[j])

    # sample time constraints
    m.initial_time_rule = pyo.Constraint(rule=initial_time)
    m.time_transition_rule = pyo.Constraint(rule=time_transition)

    # env condition constraints
    m.initial_env_rule = pyo.Constraint(range(0,2), rule=initial_env)
    m.calc_env_cond_rule = pyo.Constraint(range(0,2), rule=calc_env_cond)
    if obj == 'track_demand':
        m.calc_reference_rule = pyo.Constraint(range(0,2), rule=calc_reference)

    # physics constraints
    m.initial_state_rule = pyo.Constraint(range(0,2), rule=initial_state)
    m.state_transition_rule = pyo.Constraint(range(0,2), rule=state_transition)

    # ouput and input bounds
    m.output1_bound_rule_lower = pyo.Constraint(range(0,len(m.nn1.outputs)), rule=output1_lower_bound)
    m.output2_bound_rule_lower = pyo.Constraint(range(0,len(m.nn2.outputs)), rule=output2_lower_bound)
    m.inputs1_rule_lower = pyo.Constraint(range(0,len(m.nn1.inputs)), rule=input1_rule_lower)
    m.inputs2_rule_lower = pyo.Constraint(range(0,len(m.nn2.inputs)), rule=input2_rule_lower)
    m.inputs1_rule_upper = pyo.Constraint(range(0,len(m.nn1.inputs)), rule=input1_rule_upper)
    m.inputs2_rule_upper = pyo.Constraint(range(0,len(m.nn2.inputs)), rule=input2_rule_upper)

    # define the objective function
    if obj == 'combined': # J = demand + (-1 x productivity)
        w_p = 1.0; w_d = 1.0
        m.obj = pyo.Objective(expr=((-1*w_p*m.nn1.outputs[2])+(w_d*m.nn1.outputs[3]) +\
                                    (-1*w_p*m.nn2.outputs[2])+(w_d*m.nn2.outputs[3])) / Nh)
    if obj == 'min_demand':  # J = demand 
        m.obj = pyo.Objective(expr=(m.nn1.outputs[3] + m.nn2.outputs[3]) / Nh)
    if obj == 'max_productivity':  # J = -1 x productivity
        m.obj = pyo.Objective(expr=((-1*m.nn1.outputs[2]) + (-1*m.nn2.outputs[2])) / Nh)
    if obj == 'track_demand':  # J = (demand - reference)^2
        m.obj = pyo.Objective(expr=((m.nn1.outputs[3] - m.ref[0])**2 + (m.nn2.outputs[3] - m.ref[1])**2) / Nh)
    if obj == 'minimize_net':  # J = (weight x demand) + (intensity x demand - productivity)
        w_p = 0.0; w_d = 1.0
        m.obj = pyo.Objective(expr=(w_d*(m.nn1.outputs[3]) + w_p*(m.intensity[0]*m.nn1.outputs[3]-m.nn1.outputs[2]) +\
                                    w_d*(m.nn2.outputs[3]) + w_p*(m.intensity[1]*m.nn2.outputs[3]-m.nn2.outputs[2])) / Nh)
    if obj == 'minimize_cost':  # J = (cost x demand) + (intensity x demand - productivity)
        m.obj = pyo.Objective(expr=((m.cost[0]*m.nn1.outputs[3]) + (m.intensity[0]*m.nn1.outputs[3]-m.nn1.outputs[2]) +\
                                    (m.cost[1]*m.nn2.outputs[3]) + (m.intensity[1]*m.nn2.outputs[3]-m.nn2.outputs[2])) / Nh)
    
    return m

def test_model(dac, m, filename_test, scaler, Nh, obj, ref_signal):

    dac.eval()

    # load the real environmental data for testing
    df = pd.read_csv(filename_test)
    df["Relative Humidity (%)"] = df["Relative Humidity (%)"] / 100 # convert from percentage to decimal

    lower_bounds = {"Temperature (degC)" : 5.0, "Relative Humidity (%)" : 0.1}
    upper_bounds = {"Temperature (degC)" : 40.0, "Relative Humidity (%)" : 0.9}
    lower = pd.Series(lower_bounds)
    upper = pd.Series(upper_bounds)
    df = df.clip(lower=lower, upper=upper, axis=1)

    env_test = df[["Temperature (degC)", "Relative Humidity (%)"]].to_numpy()
    env_test_norm = scaler.transform(np.hstack((env_test, np.zeros((env_test.shape[0], 11)))))
    
    if obj == 'track_demand':
        ref = build_reference(obj)
        x_norm = np.hstack((env_test_norm[:,:2], ref))
    elif obj == 'minimize_net':
        intensity = build_reference(obj)
        x_norm = np.hstack((env_test_norm[:,:2], intensity))
    elif obj == 'minimize_cost':
        intensity_cost = build_reference(obj)
        x_norm = np.hstack((env_test_norm[:,:2], intensity_cost))
    else:
        x_norm = env_test_norm[:,:2]

        # resample to 3 times a day (number of DAC cycles in a day)
        x_norm = x_norm #x_norm[::8,:]
    state_base = torch.ones(2, dtype=torch.float32)*0.0
    state_opt = torch.ones(2, dtype=torch.float32)*0.0

    trajectories = {
        "ref_all_base" : [],
        "ref_all_opt" : [],
        "env_1_all" : [],
        "env_2_all" : [],
        "ctrl_1_all_base" : [],
        "ctrl_2_all_base" : [],
        "ctrl_3_all_base" : [],
        "ctrl_4_all_base" : [],
        "ctrl_5_all_base" : [],
        "ctrl_6_all_base" : [],
        "state1_all_base" : [],
        "state2_all_base" : [],
        "productivity_all_base" : [],
        "demand_all_base" : [],
        "t_cycle_cum_base" : [],
        "ctrl_1_all_opt" : [],
        "ctrl_2_all_opt" : [],
        "ctrl_3_all_opt" : [],
        "ctrl_4_all_opt" : [],
        "ctrl_5_all_opt" : [],
        "ctrl_6_all_opt" : [],
        "state1_all_opt" : [],
        "state2_all_opt" : [],
        "productivity_all_opt" : [],
        "demand_all_opt" : [],
        "t_cycle_cum_opt" : []
    }

    solver = pyo.SolverFactory('ipopt')
    # solver.options['max_iter'] = 100
    
    Nh_reg = 12
    t_curr = 0
    t_curr_base = 0
    t_curr_opt = 0
    t1_tot = tm.time()
    # for s, x_var in enumerate(x_norm):
        # print(f"Sample {s}/{x_norm.shape[0]}", end="\r")
    while np.min([t_curr, t_curr_base, t_curr_opt]) < x_norm.shape[0]-Nh_reg-1:
        print(f"time (abs {t_curr}, base: {t_curr_base}, opt: {t_curr_opt}) ", end="\r")
        
        if t_curr_base >= 8760:
            t_curr_base = 8759
        if t_curr_opt >= 8760:
            t_curr_opt = 8759
        x_var_base = x_norm[t_curr_base]
        x_var_opt = x_norm[t_curr_opt]

        # fix the reference signal variable (if exists)
        if (obj == 'track_demand') or (obj == 'minimize_cost') or (obj == 'minimize_net'):            
            env_var_base = x_var_base[:2]
            env_var_opt = x_var_opt[:2]

            ref_var_base = x_var_base[-1]
            ref_var_opt = x_var_opt[-1]
        else:
            env_var_base = x_var_base
            env_var_opt = x_var_opt

        # baseline control
        if t_curr_base < x_norm.shape[0]-Nh_reg-1:
            with torch.no_grad():
                u_star_norm_base = np.array([0.028,0.05,100,2400,15000,8000])
                param_vec_base = np.array([1.43,0.1,1144])
                in_vec = np.concatenate((np.zeros_like(env_var_base), np.zeros_like(state_base), u_star_norm_base, param_vec_base))
                in_vec_norm_tmp = scaler.transform(in_vec.reshape(1,-1)).flatten()
                u_star_norm_base_norm = in_vec_norm_tmp[4:10]
                param_vec_base_norm = in_vec_norm_tmp[10:]
                in_vec_norm = np.concatenate((env_var_base, state_base, u_star_norm_base_norm, param_vec_base_norm))
                in_vec_torch = torch.tensor(in_vec_norm, dtype=torch.float32)

                out_star_base = dac(in_vec_torch)
            
            # record values for plotting
            trajectories["ctrl_1_all_base"].append(in_vec[4])
            trajectories["ctrl_2_all_base"].append(in_vec[5])
            trajectories["ctrl_3_all_base"].append(in_vec[6])
            trajectories["ctrl_4_all_base"].append(in_vec[7])
            trajectories["ctrl_5_all_base"].append(in_vec[8])
            trajectories["ctrl_6_all_base"].append(in_vec[9])
            trajectories["state1_all_base"].append(out_star_base.numpy()[0])
            trajectories["state2_all_base"].append(out_star_base.numpy()[1])
            trajectories["productivity_all_base"].append(out_star_base.numpy()[2])
            trajectories["demand_all_base"].append(out_star_base.numpy()[3])
            # update the time step based on the cycle time
            trajectories["t_cycle_cum_base"].append(t_curr_base)
            t_curr_base = t_curr_base + int(np.ceil(out_star_base[4]))
            # update the states for the next cycle
            state_base = torch.tensor(out_star_base.numpy()[0:2])

        ### SOLVE OPTIMIZATION WITH NEW PARAMETERS ###
        if t_curr_opt < x_norm.shape[0]-Nh_reg-1:
            #### --- perform regression to get continuous 2nd order polynomial --- ###
            x = np.arange(0,Nh_reg,1).reshape(-1,1)
            poly = PolynomialFeatures(degree=2, include_bias=True) # Set include_bias=False if LinearRegression's fit_intercept is True (default)
            xcol = poly.fit_transform(x)
            lsr_1 = LinearRegression()
            env_1 = x_norm[t_curr_opt:t_curr_opt+Nh_reg,0].reshape(-1,1)
            lsr_1.fit(xcol, env_1)
            lsr_2 = LinearRegression()
            env_2 = x_norm[t_curr_opt:t_curr_opt+Nh_reg,1].reshape(-1,1)
            lsr_2.fit(xcol, env_2)
            beta3 = torch.tensor([lsr_1.intercept_[0],lsr_2.intercept_[0]], dtype=torch.float32)
            beta2 = torch.tensor([lsr_1.coef_[0,1],lsr_2.coef_[0,1]], dtype=torch.float32)
            beta1 = torch.tensor([lsr_1.coef_[0,2],lsr_2.coef_[0,2]], dtype=torch.float32)
            [m.beta3[i].fix(float(beta3[i])) for i in range(len(beta3))]
            [m.beta2[i].fix(float(beta2[i])) for i in range(len(beta2))]
            [m.beta1[i].fix(float(beta1[i])) for i in range(len(beta1))]
            if obj == 'track_demand':
                ref_tmp = x_norm[t_curr_opt:t_curr_opt+Nh_reg,2].reshape(-1,1)
                Ts_0 = 0
                K = ref_tmp.max()
                C = ref_tmp.min()
                k = -10
                for j in range(Nh_reg-1):
                    if (ref_tmp[j] - ref_tmp[j+1]) > 0:
                        Ts_0 = j+0.5
                        k = 10
                    if (ref_tmp[j] - ref_tmp[j+1]) < 0:
                        Ts_0 = j+0.5
                        k = -10
                m.Ts_0.fix(float(Ts_0))
                m.K.fix(float(K))
                m.C.fix(float(C))
                m.k.fix(float(k))
            if obj == 'minimize_net':
                intensity_tmp = x_norm[t_curr_opt:t_curr_opt+Nh_reg,2].reshape(-1,1)
                [m.intensity[i].fix(intensity_tmp[j]) for i, j in enumerate(np.array([0,-1]))]
            if obj == 'minimize_cost':
                Nh_c = 6
                intensity_tmp = x_norm[t_curr_opt:t_curr_opt+Nh_c,2].reshape(-1,1)
                [m.intensity[i].fix(intensity_tmp[j]) for i, j in enumerate(np.array([0,-1]))]
                cost_tmp = x_norm[t_curr_opt:t_curr_opt+Nh_c,3].reshape(-1,1)
                [m.cost[i].fix(cost_tmp[j]) for i, j in enumerate(np.array([0,-1]))]

            # fix the environmental variables
            [m.env_vars0[i].fix(float(env_var_opt[i])) for i in range(len(env_var_opt))]
            # fix the initial state variables
            [m.state_0_vars[i].fix(float(state_opt[i])) for i in range(len(state_opt))]
            # solve optimization problem
            results = solver.solve(m, tee=False)
            u_star_norm_final = np.array([m.ctrl_vars[0,i].value for i in range(6)])
            out_star_opt = np.array([m.nn1.outputs[i].value for i in m.nn1.outputs._index_set])
            # un-normalize optimal values optimization problem
            input_vec_opt = scaler.inverse_transform(np.hstack((env_var_opt, state_opt.numpy(), u_star_norm_final, np.array([0.5,0.5,0.5]))).reshape(1,-1)).flatten()
            # record values for plotting
            trajectories["ctrl_1_all_opt"].append(input_vec_opt[4])
            trajectories["ctrl_2_all_opt"].append(input_vec_opt[5])
            trajectories["ctrl_3_all_opt"].append(input_vec_opt[6])
            trajectories["ctrl_4_all_opt"].append(input_vec_opt[7])
            trajectories["ctrl_5_all_opt"].append(input_vec_opt[8])
            trajectories["ctrl_6_all_opt"].append(input_vec_opt[9])
            trajectories["state1_all_opt"].append(out_star_opt[0])
            trajectories["state2_all_opt"].append(out_star_opt[1])
            trajectories["productivity_all_opt"].append(out_star_opt[2])
            trajectories["demand_all_opt"].append(out_star_opt[3])
            # update the time step based on the cycle time
            trajectories["t_cycle_cum_opt"].append(t_curr_opt)
            t_curr_opt = t_curr_opt + int(np.ceil(out_star_opt[4]))
            # update the states for the next cycle
            state_opt = torch.tensor(out_star_opt[0:2], dtype=torch.float32)

        x_var = x_norm[t_curr]
        env_var = x_var[:2]
        input_vec = scaler.inverse_transform(np.hstack((env_var, np.zeros(11))).reshape(1,-1)).flatten()
        trajectories["env_1_all"].append(input_vec[0])
        trajectories["env_2_all"].append(input_vec[1])
        t_curr = t_curr + 1

        if (obj == 'track_demand') or (obj == 'minimize_cost') or (obj == 'minimize_net'):
            trajectories["ref_all_base"].append(float(ref_var_base))
            trajectories["ref_all_opt"].append(float(ref_var_opt))
        else:
            trajectories["ref_all_base"].append(None)
            trajectories["ref_all_opt"].append(None)

    return trajectories

def vis_test(PATH, trajectories, show=True, save=True):
    
    # plot the environmental variables
    fig1, ax = plt.subplots(2,1,sharex=True,layout="compressed")
    ax[0].plot(trajectories["env_1_all"])
    ax[0].set_ylabel("Temperature [C]")

    ax[1].plot(trajectories["env_2_all"])
    ax[1].set_ylabel("Relative Humidity")
    ax[1].set_xlabel("Time [hr]")

    # plot the optimal control variables
    num_ctrl = 6
    ctrl_names = ["u_feed_m_s","p_vac_bar","T_regenC","t_heat_s","t_des_s", "t_ads_s"]
    fig2, ax = plt.subplots(num_ctrl,1,sharex=True,layout="compressed")
    for i in range(num_ctrl):
        ax[i].plot(trajectories[f"ctrl_{i+1}_all_base"], color="b", label="base")
        ax[i].plot(trajectories[f"ctrl_{i+1}_all_opt"], color="g", label="opt")
        ax[i].set_ylabel(ctrl_names[i])
    ax[-1].set_xlabel("Time [hr]")
    ax[0].legend()

    fig3, ax = plt.subplots(4,1,sharex=True,layout="compressed")
    ax[0].plot(trajectories["t_cycle_cum_base"], trajectories[f"productivity_all_base"], color="b", label="base")
    ax[0].plot(trajectories["t_cycle_cum_opt"], trajectories[f"productivity_all_opt"], color="g", label="opt")
    ax[0].set_ylabel("Productivity")
    ax[0].legend()
    ax[1].plot(trajectories["t_cycle_cum_base"], trajectories[f"demand_all_base"], ls="--", color="b")
    ax[1].plot(trajectories["t_cycle_cum_opt"], trajectories[f"demand_all_opt"], color="g")
    ax[1].set_ylabel("Demand")
    ax[2].plot(trajectories["t_cycle_cum_base"], trajectories[f"state1_all_base"], ls="--", color="b")
    ax[2].plot(trajectories["t_cycle_cum_opt"], trajectories[f"state1_all_opt"], color="g")
    ax[2].set_ylabel("State 1")
    ax[3].plot(trajectories["t_cycle_cum_base"], trajectories[f"state2_all_base"], ls="--", color="b")
    ax[3].plot(trajectories["t_cycle_cum_opt"], trajectories[f"state2_all_opt"], color="g")
    ax[3].set_ylabel("State 2")
    ax[3].set_xlabel("Time [hr]")

    # print("Control action acuracy:")
    # for i in range(num_ctrl):
    #     print(f'    {ctrl_names[i]}: {np.mean((np.array(trajectories[f"ctrl_{i+1}_all_opt"]) - np.array(trajectories[f"ctrl_{i+1}_all_pnet"]))**2)}')

    print(f"Total cycles (opt ctrl/base ctrl): {len(trajectories['demand_all_opt'])}/{len(trajectories['demand_all_base'])}")
    print(f"Total demamnd (opt ctrl/base ctrl): {np.sum(trajectories['demand_all_opt'])/len(trajectories['demand_all_opt']):.3f}/{np.sum(trajectories['demand_all_base'])/len(trajectories['demand_all_base']):.3f}")
    print(f"Total productivity (opt ctrl/base ctrl): {np.sum(trajectories['productivity_all_opt'])/len(trajectories['demand_all_opt']):.3f}/{np.sum(trajectories['productivity_all_base'])/len(trajectories['demand_all_base']):.3f}")

    if show:
        plt.show()
    if save:
        df = pd.DataFrame.from_dict(trajectories, orient='index').T
        df.to_csv(os.path.join(PATH, 'paper', 'results', 'New_DAC_MC_ctrltest_hzn_WY.csv')) #'New_DAC_MD_test.csv','New_DAC_MP_test.csv','New_DAC_CO_test.csv','New_DAC_TD_test.csv'