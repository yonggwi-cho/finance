#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import csv
import BSM as bsm

parser = argparse.ArgumentParser()
parser.add_argument("-Ns","--Nsample",type=int,default=10)
parser.add_argument("-m","--mu",type=float,required=True)
parser.add_argument("-s","--sigma",type=float,required=True)
parser.add_argument("-dt","--dt",type=float,default=1.0)
parser.add_argument("-c","--csvfile",required=True)
args = parser.parse_args()

z=list()
if (args.Nsample % 2 != 0):
    print("Nsample must be even number.")
    sys.exit(1)

# read csv file
with open(args.csvfile,"r") as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        z.append(float(row[4]))

Nt = len(z)
init = z[0]

# calc BSM model
fig = plt.figure(figsize=(16,10))
t=range(Nt)
for i in range(args.Nsample):
    model = bsm.BSM(Nt,args.mu,args.sigma,t_init=0.0,y_init=init)
    y = model.predict()
    rows = args.Nsample/2
    cols = args.Nsample/rows
    plt.subplot(cols,rows,i+1)
    # plot BTC data
    plt.plot(t,z,"bs:",markeredgecolor="black")
    # plot BSM model
    plt.plot(t,y,"ro-",markeredgecolor="black")
    plt.xlabel("time")
    plt.ylabel("BTC")
    plt.yscale("log")
    del model
fig.suptitle("(mu,sigma)=("+str(args.mu)+","+str(args.sigma)+")",fontsize=30)
plt.show()
