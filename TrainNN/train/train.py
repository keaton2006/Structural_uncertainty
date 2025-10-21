# Uncertainty in Physics-Informed Inverse Problems
# Anonymous, Anonymous 5/2025
# This code was originally developed
# from the code produced by the following project.
# Hamiltonian Neural Networks | 2019
# Sam Greydanus, Misko Dzamba, Jason Yosinski

import torch, argparse
import numpy as np
from sklearn.model_selection import train_test_split

import os, sys
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PARENT_DIR)

from nn_models import MLP
from hnn import HNN
from data import get_dataset
from utils import L2_loss
import sys


def get_args():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('--time_index', default=0, type=int, help='time index')
    parser.add_argument('--input_dim', default=2, type=int, help='dimensionality of input tensor')
    parser.add_argument('--hidden_dim', default=100, type=int, help='hidden dimension of mlp')
    parser.add_argument('--learn_rate', default=1e-2, type=float, help='learning rate')
    parser.add_argument('--nonlinearity', default='tanh', type=str, help='neural net nonlinearity')
    parser.add_argument('--total_steps', default=400000, type=int, help='number of gradient steps')
    parser.add_argument('--print_every', default=1, type=int, help='number of gradient steps between prints')
    parser.add_argument('--name', default='training-results', type=str, help='only one option right now')
    parser.add_argument('--baseline', dest='baseline', action='store_true', help='run baseline or experiment?')
    parser.add_argument('--use_rk4', dest='use_rk4', action='store_true', help='integrate derivative with RK4')
    parser.add_argument('--verbose', dest='verbose', action='store_true', help='verbose?')
    parser.add_argument('--field_type', default='solenoidal', type=str, help='type of vector field to learn')
    parser.add_argument('--seed', default=8, type=int, help='random seed')
    parser.set_defaults(feature=True)
    return parser.parse_args()

def train(args, time_index):
  print(time_index)
  # set random seed
  torch.manual_seed(args.seed)
  np.random.seed(args.seed)

  # init model and optimizer
  if args.verbose:
    print("Training baseline model:" if args.baseline else "Training HNN model:")

  output_dim = args.input_dim if args.baseline else 2
  nn_model = MLP(args.input_dim, args.hidden_dim, output_dim, args.nonlinearity)
  model = HNN(args.input_dim, differentiable_model=nn_model,
              field_type=args.field_type, baseline=args.baseline)
  optim = torch.optim.Adam(model.parameters(), args.learn_rate, weight_decay=1e-4)

  data = np.load("../traindata/TrainHNN_GZ_26_origin_1.npy")
  data = data[time_index][150:250,100:200]
  data = data.reshape([100*100,7])
  print()
  data_x, data_test_x, data_dIdt, data_test_dIdt = train_test_split(data[:,:7], data[:,3], test_size=0.1, random_state=0)
  x = torch.tensor(data_x, requires_grad=True, dtype=torch.float32)
  
  test_x = torch.tensor( data_test_x, requires_grad=True, dtype=torch.float32)
  dIdt = torch.Tensor(data_dIdt)
  test_dIdt = torch.Tensor(data_test_dIdt)

  # vanilla train loop
  stats = {'train_loss': [], 'test_loss': []}
  print(args.verbose)
  for step in range(args.total_steps+1):
    if step==100000:
      optim.param_groups[0]['lr'] = default=1e-3
      pg = optim.param_groups[0]
      print({i: pg[i] for i in list(pg.keys())[1:]})
    elif step==200000:
      optim.param_groups[0]['lr'] = default=1e-4
      pg = optim.param_groups[0]
      print({i: pg[i] for i in list(pg.keys())[1:]})
    elif step==300000:
      optim.param_groups[0]['lr'] = default=1e-5
      pg = optim.param_groups[0]
      print({i: pg[i] for i in list(pg.keys())[1:]})
    
    # train step
    dxdt_hat = model.rk4_time_derivative(x[:,1:3]) if args.use_rk4 else model.time_derivative(x[:,1:3])
    x_sym = data_x[:,1:3]
    x_sym[:,1] = -x_sym[:,1]
    x_sym = torch.tensor(x_sym, requires_grad=True, dtype=torch.float32)
    h_hat1, h_hat = model.forward(x[:,1:3])
    h_hat_sym1, h_hat_sym = model.forward(x_sym)

    loss_sym = 10000*L2_loss(h_hat.reshape(9000), h_hat_sym.reshape(9000))
    delta_k = 3.0
    D = 3.0
    delta_omega = 0.4
    ky = 1.0
    gamma = (ky*(x[:,2]**2+ky**2))/(D*(1+x[:,2]**2+ky**2)**3)*torch.exp(-(x[:,2]/delta_k)**2)
    C = gamma*x[:,6] - delta_omega*x[:,6]**2
    dIdt_hat = -dxdt_hat[:,0]*x[:,4] + dxdt_hat[:,1]*x[:,5] + C
    loss_l2 = L2_loss(10000*dIdt, 10000*dIdt_hat) 

    loss = loss_l2 + loss_sym
    
    loss.backward() ; optim.step() ; optim.zero_grad()
    
    test_dxdt_hat = model.rk4_time_derivative(test_x[:,1:3]) if args.use_rk4 else model.time_derivative(test_x[:,1:3])

    gamma = (ky*(test_x[:,2]**2+ky**2))/(D*(1+test_x[:,2]**2+ky**2)**3)*torch.exp(-(test_x[:,2]/delta_k)**2)
    C = gamma*test_x[:,6] - delta_omega*test_x[:,6]**2
    test_dIdt_hat = -test_dxdt_hat[:,0]*test_x[:,4] + test_dxdt_hat[:,1]*test_x[:,5] + C
    test_loss = L2_loss(test_dIdt, test_dIdt_hat)

    # logging
    stats['train_loss'].append(loss.item())
    stats['test_loss'].append(test_loss.item())
    if args.verbose and step % args.print_every == 0:
      print("step {}, train_loss {:.10e}, test_loss {:.10e}".format(step, loss.item(), test_loss.item()))

  return model, stats

if __name__ == "__main__":
    args = get_args()
    args.verbose = True
    time_index = args.time_index
    model, stats = train(args,time_index)

    # save
    os.makedirs(args.save_dir) if not os.path.exists(args.save_dir) else None
    path = '../results/{}_{}.tar'.format(args.name, time_index)
    torch.save(model.state_dict(), path)
