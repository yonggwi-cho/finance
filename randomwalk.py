#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

argvs=sys.argv
N=int(argvs[1])
Init=int(argvs[2])
mu=float(argvs[3])
sigma=float(argvs[4])

x=np.array([ i for i in range(N)])
b=np.zeros(N)
y=np.zeros(N)

rand = np.random.normal(size=(N))
for it in range(N):
    if it == 1:
        b[it] = 1
    else :
        b[it] = b[it-1] + rand[it]
    y[it] = Init + mu*it + sigma*b[it]
        
plt.plot(x,y,"b")
plt.show()
