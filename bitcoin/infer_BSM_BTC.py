#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import csv
import optimizer as opt


if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    #parser.add_argument("-N","--N",default=100)
    parser.add_argument("-Ns","--Nsample",default=100)
    #parser.add_argument("-I","--init",default=1.0)
    parser.add_argument("-m","--mu",type=float,default=0.4)
    parser.add_argument("-s","--sigma",type=float,default=0.1)
    parser.add_argument("-seed","--seed",type=int,default=1)
    parser.add_argument("-c","--csvfile",type=str,default="./BTCJPY_test.csv")
    args = parser.parse_args()

    # calc loss for BSM model and BTC data
    optimizer = opt.Optimizer(args.Nsample,args.csvfile)
    ave, var = optimizer.estimate_model(args)
    print("averaged_squared_error= "+str(ave))
    print("RMSerror= "+str(math.sqrt(ave)))
    print("variance= "+str(var))
