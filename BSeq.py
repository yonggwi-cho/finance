#!/usr/bin/env python
# diffusion equation on 1D space
# Black-Showles equation

import numpy as np
import pylab as pl
import math
from mpl_toolkits.mplot3d import Axes3D
import sys

# params
Nt=150
Ns=200
Tmax=1
r= 0.1
sgm=0.3
K=100
alpha=0.5-r/sgm**2
beta=-(sgm**2+2*r)**2/(8*sgm**2)
xmin=-5
xmax=5
dx=float(xmax-xmin)/Ns
dtau=0.5*sgm**2/Tmax/Nt
d=dtau/(dx)**2

# initialize
U = np.zeros((Ns,Nt))
X = np.zeros(Ns)
Tau = np.zeros(Nt)
for ix in range(Ns):
    X[ix] = -5 + dx*ix
for itau in range(Nt):
    Tau[itau] = 0 + dtau*itau
print Tau

# I.C
for i in range(Ns):
    if X[i] >= math.log(K) :
        U[i][0] = math.exp((1-alpha)*X[i])-K*math.exp(-alpha*X[i])
    else:
        U[i][0] = 0.0

# B.C
for j in range(Nt):
    U[0][j] = 0.0
    U[Ns-1][j] = math.exp(X[Ns-1])*math.exp(-alpha*X[Ns-1]-beta*Tau[j])

# time evolution
for j in range(Nt-1):
    for i in range(Ns-1):
        U[i][j+1] = d*U[i+1][j] +(1.0-2.0*d)*U[i][j] +d*U[i-1][j]

# write data
output = open("./output.dat","w")
for i in range(Ns):
#    for j in range(Nt):
    for j in [Nt-1]:
        # tau -> t
        t=Tmax-2.0*Tau[j]/sgm**2
        # x -> S
        S = math.exp(X[i])
        # U -> V
        V = math.exp(alpha*X[i]+beta*Tau[j])*U[i][j]
        output.write(str(S)+"\t"+str(t)+"\t"+str(V)+"\n")
output.close()
