#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-N","--N",default=100)
parser.add_argument("-I","--init",default=1.0)
parser.add_argument("-m","--mu",type=float,default=1.0)
parser.add_argument("-s","--sigma",type=float,default=1.0)
parser.add_argument("-dt","--dt",type=float,default=1.0)
args = parser.parse_args()

x=np.array([ i for i in range(args.N)])
w=np.zeros(args.N,dtype=float)
y=np.zeros(args.N,dtype=float)
z=np.zeros(args.N,dtype=float)
t=np.zeros(args.N,dtype=float)

y0 = float(args.init)
t0 = 0
t[0] = 0
y[0] = args.init


for it in range(1,args.N):
    t += args.dt
    if it == 0:
        #w[it] = np.random.normal()
        w[it] = 0.0
    else :
        w[it] = w[it-1] + np.random.normal()

    t[it] = t[it-1] + args.dt
    y[it] = y[0]*math.exp((args.mu-0.5*args.sigma**2)*t[it]+args.sigma*w[it])

    print y[it]

plt.plot(t,y,"b")
plt.yscale("log")
plt.show()
