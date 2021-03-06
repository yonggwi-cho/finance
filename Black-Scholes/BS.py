#!/usr/bin/env python
# BS equation 

import sys
import math
import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D

# init condition
def CallPayoff(x, k):
    y = math.exp(0.5*(k+1)*x)-math.exp(0.5*(k-1)*x) # ?
    if y < 0 :
        #print "y=",y
        y=0.0
    return y

# Boundary condition at S=0 
def CallBC1(x,tau):
    return 0.0

# Boundary condition for S -> infinity
def CallBC2(x,tau,alpha,beta,k):
    return math.exp((1.0-alpha)*x-beta*tau)-math.exp(-(k+beta)*tau)

# define paramaters
argvs=sys.argv
r=float(argvs[1])
sigma=float(argvs[2])
N=200
M=150
xmax=5.0
T=1.0
K=100
k=2.0*r/sigma**2
alpha=-1.0/2.0*(k-1)
beta=-1.0/4.0*(k+1)**2 # ?
dt = T*(1.0/2.0)*sigma**2/M
dx = xmax/N
a = dt/dx**2

# prepare numpy array
U0 = np.zeros((2*N+1),float)
U  = np.zeros((2*N+1),float)
V  = np.zeros((2*N+1),float)
X  = np.zeros((2*N+1),float)

# initial values
for j in range(-N,N+1):
    U0[j+N] = CallPayoff(j*dx,k)

for i in range(M):
    tau = i*dt
    # Boundary condition
    U[0] = CallBC1(-N*dx,tau)
    U[2*N] = CallBC2(N*dx,tau,alpha,beta,k)

    # time evolution
    for j in range(-N+1,N):
        U[j+N] = U0[j+N] + a*(U0[j+N-1]-2*U0[j+N]+U0[j+N+1])

    for j in range(-N,N+1):
        U0[j+N] = U[j+N]

for j in range(-N,N+1):
    V[j+N] = K*math.exp(alpha*j*dx+beta*0.5*sigma**2*T)*U0[j+N]

# write data
output = open("./output_BS.dat","w")
output.write("dt/dx**2 = "+str(a))
for j in range(2*N+1):
        output.write(str(K*math.exp(xmax/N*(j-N)))+"\t"+str(V[j])+"\n")
        X[j] = K*math.exp(xmax/N*(j-N))
output.close()

# plot 
plt.plot(X,V,"r-")
plt.show()
