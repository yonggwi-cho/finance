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
parser.add_argument("-m","--mu",type=float,default=1.0)
parser.add_argument("-s","--sigma",type=float,default=1.0)
parser.add_argument("-dt","--dt",type=float,default=1.0)
parser.add_argument("-c","--csvfile",required=True)
parser.add_argument("-se","--seed",type=int,default=-1)
args = parser.parse_args()

z=list()

# read csv file
with open(args.csvfile,"r") as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        z.append(float(row[4]))

Nt = len(z)
init = z[0]

# calc BSM model
model = bsm.BSM(Nt,args.mu,args.sigma,y_init=init,seed=args.seed)
y = model.predict()

# plot
t = range(Nt)
# plot BSM model
plt.plot(t,y,"ro-",label="BSM model",markeredgecolor="black")
# plot BTC data
plt.plot(t,z,"bs:",label="market price",markeredgecolor="black")
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("BTC")
plt.legend()
plt.show()
