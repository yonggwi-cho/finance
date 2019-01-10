#!/usr/bin/env python
# diffusion equation on 1D space

import numpy as np
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

N=200
M=150
xmax=5.0
r=0.3


T = np.zeros((Ns,Nt))

# initial values
for i in range(Ns):
    T[i][0] = 100.0

# Boundary condition
for j in range(Nt):
    T[0][j] = 0.0
    T[Ns-1][j] = 0.0

# time evolution
for j in range(Nt-1):
    for i in range(Ns-1):
        T[i][j+1] = T[i][j] + 0.5*(T[i+1][j]+T[i-1][j]-2.0*T[i][j])

# write data
output = open("./output.dat","w")
for i in range(Ns):
    for j in range(Nt):
        output.write(str(i)+"\t"+str(j)+"\t"+str(T[i][j])+"\n")
output.close()

