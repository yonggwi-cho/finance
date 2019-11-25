#/!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import commands as cm

# sampling num
Nsample = 10

# initial values
mu = 0.01
sgm = 0.1

# list
loss = list()
var  = list()
mu_list = list()
sgm_list = list()

# gaussian sampler for mu
for i in range(1,Nsample+1):
    mup = i*0.1*mu
    output=cm.getoutput("python train_BSM_BTC.py -m "+str(mup)+" -s "+str(sgm)).split("\n")
    loss = output[0].split(" ")[1]
    var = output[1].split(" ")[1]
    print(str(mup)+"\t"+str(sgm)+"\t"+str(loss)+"\t"+str(var))

# gaussian sampler for sgm
for i in range(1,Nsample+1):
    sgmp = i*0.1*sgm
    output=cm.getoutput("python train_BSM_BTC.py -m "+str(mu)+" -s "+str(sgmp)).split("\n")
    loss = output[0].split(" ")[1]
    var = output[1].split(" ")[1]
    print(str(mu)+"\t"+str(sgmp)+"\t"+str(loss)+"\t"+str(var))

# gaussian sampler for mu and sgm
for j in range(1,Nsample+1):
    for i in range(1,Nsample+1):
        mup = i*0.1*mu
        sgmp = j*0.1*sgm
        output=cm.getoutput("python train_BSM_BTC.py -m "+str(mup)+" -s "+str(sgmp)).split("\n")
        loss = output[0].split(" ")[1]
        var = output[1].split(" ")[1]
        print(str(mup)+"\t"+str(sgmp)+"\t"+str(loss)+"\t"+str(var))
