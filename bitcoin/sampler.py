#/!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import commands as cm

# sampling num
Nsample = 100

# initial values
mu = 0.5
sgm = 0.1

# list
loss = list()
var  = list()
mu_list = list()
sgm_list = list()

# gaussian sampler for mu
for i in range(Nsample):
    mup = np.random.normal(mu,mu/10.0)
    output=cm.getoutput("python train_BSM_BTC.py -m "+str(mup)+" -s "+str(sgm)).split("\n")
    loss = output[0].split(" ")[1]
    var = output[1].split(" ")[1]
    print(str(mup)+"\t"+str(sgm)+"\t"+str(loss)+"\t"+str(var))

# gaussian sampler for sgm
for i in range(Nsample):
    sgmp = np.random.normal(sgm,sgm/10.0)
    output=cm.getoutput("python train_BSM_BTC.py -m "+str(mu)+" -s "+str(sgmp)).split("\n")
    loss = output[0].split(" ")[1]
    var = output[1].split(" ")[1]
    print(str(mu)+"\t"+str(sgmp)+"\t"+str(loss)+"\t"+str(var))

# gaussian sampler for mu and sgm
for i in range(Nsample):
    mup = np.random.normal(mu,mu/10.0)
    sgmp = np.random.normal(sgm,sgm/10.0)
    output=cm.getoutput("python train_BSM_BTC.py -m "+str(mup)+" -s "+str(sgmp)).split("\n")
    loss = output[0].split(" ")[1]
    var = output[1].split(" ")[1]
    print(str(mup)+"\t"+str(sgmp)+"\t"+str(loss)+"\t"+str(var))
