#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-N","--N",default=1000)
parser.add_argument("-I","--init",default=1.0)
parser.add_argument("-m","--mu",type=float,default=1.0)
parser.add_argument("-s","--sigma",type=float,default=1.0)
parser.add_argument("-dt","--dt",type=float,default=1.0)
args = parser.parse_args()

x=np.array([ i for i in range(args.N)])
w=np.zeros(args.N)
y=np.zeros(args.N)
z=np.zeros(args.N)

y0 = args.init
t0 = 0
t = 0
y[0] = args.init

for it in range(1,args.N):
    t += args.dt
    if it == 0:
        w[it] = np.random.normal()
    else :
        w[it] = w[it-1] + np.random.normal()

    y[it] = y[it-1] + args.mu*y[it-1]*args.dt + args.sigma*y[it-1]*w[it-1]
    print y[it]

plt.plot(x,y,"b")
plt.yscale("log")
plt.show()
