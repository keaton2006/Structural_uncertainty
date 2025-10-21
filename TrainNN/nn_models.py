# Uncertainty in Physics-Informed Inverse Problems
# Anonymous, Anonymous 2025
# This code was originally developed
# from the code produced by the following project.
# Hamiltonian Neural Networks | 2019
# Sam Greydanus, Misko Dzamba, Jason Yosinski

import torch
import numpy as np
from utils import choose_nonlinearity

class MLP(torch.nn.Module):
  def __init__(self, input_dim, hidden_dim, output_dim, nonlinearity='tanh'):
    super(MLP, self).__init__()
    self.linear1 = torch.nn.Linear(input_dim, 100)
    self.linear2 = torch.nn.Linear(100, 10)
    self.linear3 = torch.nn.Linear(10, output_dim, bias=None)

    for l in [self.linear1, self.linear2, self.linear3]:
      torch.nn.init.orthogonal_(l.weight) # use a principled initialization

    self.nonlinearity = choose_nonlinearity(nonlinearity)

  def forward(self, x, separate_fields=False):
    h = self.nonlinearity( self.linear1(x) )
    h = self.nonlinearity( self.linear2(h) )
    return self.linear3(h)
