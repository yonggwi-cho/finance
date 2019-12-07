#!/usr/bin/env python

# library
import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
import math

# variables
N_VC = 10
Nt = 1000

# mulkov process
# v[t+1] =  alpha*v[t] + beta*random[t]
v = np.zeros((Nt,N_VC))
#for i in range(N_VC):
#    v[0][i] = rnd.random(1)*10

# model parameters
alpha = np.ones((N_VC))
beta  = np.ones((N_VC))

def update_model():
    alpha = rnd.random(N_VC)
    beta  = rnd.random(N_VC)

def step():
    for it in range(Nt-1):
        for i in range(N_VC):
            v[it+1,i] = alpha[i]*v[it,i] + pow(-1.0,rnd.choice([1,2]))*beta[i]*rnd.random()

#def env():

#def train():

def plot():
    x = range(Nt)
    y = v[:,0]
    for i in range(N_VC):
        plt.plot(x,v[:,i])
    plt.pause(1)

if __name__ == "__main__" :
    #train()
    budget = 1000
    while budget >= 0.0 :
        update_model()
        step()
        plot()
        budget -= 100
    
    
