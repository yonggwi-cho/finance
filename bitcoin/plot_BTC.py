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
        y1.append(row[1])
        y2.append(row[2])
        y3.append(row[3])
        y4.append(row[4])

# plot
x=range(len(y1))
#plt.plot(x,y1,"ro-")
#plt.plot(x,y2,"bs--")
#plt.plot(x,y3,"g^:")
plt.plot(x,y4,"mo--")
plt.yscale("log")
plt.show()
