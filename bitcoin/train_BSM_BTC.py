#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import csv
import BSM as bsm

parser = argparse.ArgumentParser()
parser.add_argument("-N","--N",default=100)
#parser.add_argument("-I","--init",default=1.0)
parser.add_argument("-m","--mu",type=float,default=0.5)
parser.add_argument("-s","--sigma",type=float,default=0.1)
parser.add_argument("-seed","--seed",type=int,default=1)
parser.add_argument("-c","--csvfile",type=str,default="./BTCJPY_train.csv")
args = parser.parse_args()

timestamp=list()
z=list()

# read csv file
with open(args.csvfile,"r") as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        timestamp.append(row[0])
        z.append(row[4])

N = len(z) # use close market's value
t_init = 0
y_init = z[t_init]

# BSM model
model = bsm.BSM(N,args.mu,args.sigma,t_init,y_init,args.seed)
t, y = model.predict_fixrandom()

# calc loss function


# plot BSM model
plt.plot(t,y,"ro-")
plt.yscale("log")

# plot BTC data
x=range(N)
plt.plot(x,z,"bs:")
plt.yscale("log")
plt.show()
