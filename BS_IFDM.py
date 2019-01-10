#!/usr/bin/env python
# BS equation
# implict FDM

import sys
import math
import numpy as np
import pylab as pl

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

# solve inverse matrix
def solve_inverse(A):
    B = np.linalg.inv(A)
    # check
    print "check:residual=",float(np.ndarray.sum(A.dot(B)))-float(A.shape[0])
    return B

# solve equation
def solve(V,B,k):
    U  = np.zeros((2*N+1)) 
    S0 = np.zeros(2*N-1)
    S = np.zeros(2*N-1)
    
    # define matrix B
    for i in range(2*N-1):
        B[i][i] = 1.0 + 2.0*a
        if i != 2*N-2:
            B[i][i+1] = -a
            B[i+1][i] = -a
        if i != 0 :
            B[i][i-1] = -a
            B[i-1][i] = -a
            
    # solve inverse of B
    print B
    Binv = solve_inverse(B)
    print Binv
    
    # I.C at tau=0
    for j in range(1,2*N):
        x = (j-N)*dx
        S0[j-1] = CallPayoff(x,k)
    for j in range(2*N+1):
        x = (j-N)*dx
        U[j] = CallPayoff(x,k)

    b = np.zeros(2*N-1)
    # time evolution
    for i in range(M):
        tau = i*dt
        
        # B.C (This does not depend on time slices!)
        U[0]   = CallBC1(-N*dx,tau)
        U[2*N] = CallBC2( N*dx,tau,alpha,beta,k)
        
        # set b vector
        b[0]     = a*U[0]
        b[2*N-2] = a*U[2*N]
        
        # mat*vec multiplication
        S = np.dot(Binv,S0+b)
                
        for j in range(2*N-1):
            S0[j]  = S[j] # for next interation
            
    for j in range(2*N-1):
        U[j+1] = S0[j] 

    # transfer U->V
    for j in range(2*N+1):
        x = (j-N)*dx
        V[j] = K*math.exp(alpha*x+beta*0.5*sigma**2*T)*U[j]

def main():
    # prepare numpy array
    B  = np.zeros((2*N-1,2*N-1))

    V  = np.zeros((2*N+1)) 
    X  = np.zeros((2*N+1)) 

    # call solve
    solve(V,B,k)

    # write data
    output = open("./output_IFDM.dat","w")
    for j in range(2*N+1):
        output.write(str(K*math.exp(xmax/N*(j-N)))+"\t"+str(V[j])+"\n")
        X[j] = K*math.exp(xmax/N*(j-N))
    output.close()

    # plot 
    pl.xlim([0,200])
    pl.ylim([0,120])
    pl.title("Black-Showles equation")
    pl.ylabel("Option Price")
    pl.xlabel("$S_0$")
    pl.plot(X,V,"ro")
    pl.show()

if __name__ == "__main__":
    # define paramaters
    argvs=sys.argv
    r=float(argvs[1])
    sigma=float(argvs[2])
    T=float(argvs[3])
    N=int(argvs[4])
    M=int(argvs[5])
    xmax=5.0
    K=100
    k=2.0*r/sigma**2
    alpha=-1.0/2.0*(k-1)
    beta=-1.0/4.0*(k+1)**2 # ?
    dt = T*(1.0/2.0)*sigma**2/M
    dx = xmax/N
    a = dt/dx**2

    # main
    main()
        
