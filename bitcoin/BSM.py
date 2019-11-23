#!/usr/bin/env python
import numpy as np
import math

class BSM():

    def __init__(self,N,mu,sgm,t_init,y_init):
        self.N = N
        self.mu = mu
        self.sgm = sgm
        self.dt = 1.0
        self.t = np.zeros(self.N,dtype=float)
        self.w = np.zeros(self.N,dtype=float) # stanrd standard browninaMotion
        self.y = np.zeros(self.N,dtype=float)
        self.t[0] = float(t_init)
        self.y[0] = float(y_init)

    def predict(self):
        # predict by BSM model
        for it in range(self.N):
            # make standard browninaMotion
            if it == 0:
                self.w[it] = 0.0
            else :
                self.w[it] = self.w[it-1] + np.random.normal()
            # BSM model
            self.y[it] = self.y[0]*math.exp((self.mu-0.5*self.sgm**2)*self.t[it]+self.sgm*self.w[it])
            self.t[it] = self.t[it-1] + self.dt
        return self.t, self.y