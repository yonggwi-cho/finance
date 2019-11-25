#!/usr/bin/env python
#-*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import csv
import argparse


parser =argparse.ArgumentParser()
parser.add_argument("-f","--file",required=True)
args = parser.parse_args()

x=list()
y1=list()
y2=list()
y3=list()
y4=list()

# read csv file
with open(args.file,"r") as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        #x.append(row[0])
        y1.append(float(row[1]))
        y2.append(float(row[2]))
        y3.append(float(row[3]))
        y4.append(float(row[4]))

# plot
x=range(len(y1))
#plt.plot(x,y1,"mo:")
#plt.plot(x,y2,"bs--")
#plt.plot(x,y3,"g^:")
plt.plot(x,y4,"ro-",markeredgecolor="black")
plt.yscale("log")
plt.xlabel("time",fontsize=20)
plt.ylabel("BTC",fontsize=20)

# save plot data
with open("plot_BTC.data","w") as file:
    file.write("#time,close\n")
    for i in range(len(x)):
        file.write(str(x[i])+"\t"+str(y4[i])+"\n")

plt.show()
