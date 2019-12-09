#/!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import commands as cm
import optimizer as opt
import argparse as ae

# get args
parser = ae.ArgumentParser()
parser.add_argument("-m","--mu",type=float,required=True)
parser.add_argument("-s","--sgm",type=float,required=True)
parser.add_argument("-Ns","--Nsample",default=10)
parser.add_argument("-c","--csvfile",type=str,default="./BTCJPY_train.csv")
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
optimizer = opt.Optimizer(args.Nsample,args.csvfile)

# initial values
optimizer.add_param(dict)
optimizer.sample()
