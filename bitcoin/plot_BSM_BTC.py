#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-N","--N",default=100)
#parser.add_argument("-I","--init",default=1.0)
parser.add_argument("-m","--mu",type=float,default=1.0)
parser.add_argument("-s","--sigma",type=float,default=1.0)
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

N = len(y4)
init = y4[0]

w=np.zeros(N,dtype=float)
y=np.zeros(N,dtype=float)
t=np.zeros(N,dtype=float)

t0 = 0
t[0] = 0.0
y[0] = init
# calc BSM model
for it in range(N):
    if it == 0:
        #w[it] = np.random.normal()
        w[it] = 0.0
    else :
        w[it] = w[it-1] + np.random.normal()
    y[it] = y[0]*math.exp((args.mu-0.5*args.sigma**2)*t[it]+args.sigma*w[it])
    t[it] = t[it-1] + args.dt


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
