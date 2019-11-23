#!/usr/bin/env python
import csv
import BSM as bsm
import numpy as np
import matplotlib.pyplot as plt
class Optimizer():

    def __init__(self,N,csvfile):
        self.Nsample = N
        self.csvfile = csvfile
        self.timestamp = list()
        self.z = list()
        self.t = list()
        self.read_csv()

    def read_csv(self):
        # read csv file
        with open(self.csvfile,"r") as file:
            reader = csv.reader(file)
            header = next(reader)
            print("read csv file...")
            print(header)
            for row in reader:
                self.timestamp.append(row[0])
                self.z.append(row[4]) # use market's close value
        self.t = range(len(self.z))

    def plot_fig(self,y):
        plt.plot(self.t,y,"ro-")    # plot BSM model
        plt.plot(self.t,self.z,"bs:")    # plot BTC data
        plt.yscale("log")
        plt.show()

    def estimate(self,args):
        # initialize
        N = len(self.z) # use close market's value
        t_init = 0
        y_init = self.z[t_init]

        # BSM model
        model = bsm.BSM(N,args.mu,args.sigma,t_init,y_init,args.seed)
        #t, y = model.predict_fixrandom()
        y = model.predict()
        loss = model.calc_loss(np.array(self.z,dtype=float))
        return loss, y
