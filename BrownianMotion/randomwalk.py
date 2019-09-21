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
args = parser.parse_args()

x=np.array([ i for i in range(args.N)])
b=np.zeros(args.N)
y=np.zeros(args.N)
z=np.zeros(args.N)

y0 = args.init
rand = np.random.normal(size=(args.N))
for it in range(args.N):
    if it == 1:
        b[it] = 1
    else :
        b[it] = b[it-1] + rand[it]

    y[it] = y0 + args.mu*it + args.sigma*b[it]

plt.plot(x,y,"b")
plt.show()
