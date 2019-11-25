#!/usr/bin/env python
import csv
import BSM as bsm
import numpy as np
import matplotlib.pyplot as plt
import math

class Optimizer():

    def __init__(self,N,csvfile):
        self.Nsample = N
        self.csvfile = csvfile
        self.timestamp = list()
        self.z = list()
        self.t = list()
        self.read_csv()
        self.seed = -1 # seed must be not setted.

    def read_csv(self):
        # read csv file
        with open(self.csvfile,"r") as file:
            reader = csv.reader(file)
            header = next(reader)
            #print("read csv file...")
            #print(header)
            for row in reader:
                self.timestamp.append(row[0])
                self.z.append(row[4]) # use market's close value
        self.t = range(len(self.z))

    def plot_fig(self,y):
        plt.plot(self.t,y,"ro-")    # plot BSM model
        plt.plot(self.t,self.z,"bs:")    # plot BTC data
        plt.yscale("log")
        plt.show()

    def estimate_last(self,args):
        # initialize
        N = len(self.z) # use close market's value
        t_init = 0
        y_init = self.z[t_init]

        # BSM model
        model = bsm.BSM(N,args.mu,args.sigma,t_init,y_init,self.seed)
        y = model.predict()
        #print(float(y[-1]),float(self.z[-1]),(float(y[-1])-float(self.z[-1])))
        error = (float(y[-1]) - float(self.z[-1]))**2
        return error

    def estimate(self,args):
        # initialize
        N = len(self.z) # use close market's value
        t_init = 0
        y_init = self.z[t_init]

        # BSM model
        model = bsm.BSM(N,args.mu,args.sigma,t_init,y_init,self.seed)
        y = model.predict()
        loss = model.calc_loss(np.array(self.z,dtype=float))
        return loss, y

    def estimate_Nsample(self,args):
        loss         = 0.0
        loss_squared = 0.0
        for i in range(self.Nsample):
            tmp, y = self.estimate(args)
            loss += tmp
            loss_squared += tmp**2

        average_loss = loss/self.Nsample
        average_loss_squread = loss_squared/self.Nsample
        variance =  average_loss_squread - average_loss

        return average_loss, math.sqrt(variance)

    def estimate_model(self,args):
        # calc error of last value
        last         = 0.0
        last_squared = 0.0
        for i in range(self.Nsample):
            tmp = self.estimate_last(args)
            last += tmp
            last_squared += tmp**2
            #print(tmp)

        average_last = last/self.Nsample
        average_last_squread = last_squared/self.Nsample
        variance =  average_last_squread - average_last

        return math.sqrt(average_last), math.sqrt(variance)

    def estimate_Nsample_wplot(self,args):
        loss         = 0.0
        loss_squared = 0.0
        for i in range(self.Nsample):
            tmp, y = self.estimate(args)
            loss += tmp
            loss_squared += tmp**2

        average_loss = loss/self.Nsample
        average_loss_squread = loss_squared/self.Nsample
        variance =  average_loss_squread - average_loss

        return average_loss, math.sqrt(variance)
