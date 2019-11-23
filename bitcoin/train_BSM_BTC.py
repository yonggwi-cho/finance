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
parser.add_argument("-dt","--dt",type=float,default=1.0)
parser.add_argument("-c","--csvfile",required=True)
args = parser.parse_args()

x=list()
y1=list()
y2=list()
y3=list()
y4=list()

# read csv file
with open(args.csvfile,"r") as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        x.append(row[0])
        y1.append(row[1])
        y2.append(row[2])
        y3.append(row[3])
        y4.append(row[4])

N = len(y4) # use close market's value
t_init = 0
y_init = y4[t_init]

#BSM model
model = bsm.BSM(N,args.mu,args.sigma,t_init,y_init)
t, y = model.predict()

# plot BSM model
plt.plot(t,y,"ro-")
plt.yscale("log")

# plot BTC data
x=range(N)
#plt.plot(x,y1,"ro-")
#plt.plot(x,y2,"bs--")
#plt.plot(x,y3,"g^:")
plt.plot(x,y4,"bs:")
plt.yscale("log")
plt.show()
