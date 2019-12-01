#/!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import commands as cm
import optimizer as opt
import argparse as ae

# get args
parser = ae.ArgumentParser()
parser.add_argument("-m","--mu",required=True)
parser.add_argument("-s","--sgm",required=True)
parser.add_argument("-Ns","--Nsample",default=10)
args = parser.parse_args()

# list
loss = list()
var  = list()
mu_list = list()
sgm_list = list()

# make dict
dict = {}
dict["mu"] = args.mu
dict["sgm"] = args.sgm
optimizer = opt.Optimizer(args.Nsample,)

# initial values
optimizer.add_param(dict[0])
optimizer.add_param(dict[1])
optimizer.sample()
