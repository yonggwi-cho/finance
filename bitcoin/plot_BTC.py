#!/usr/bin/env python
#-*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import csv
import argparse


parser =argparse.ArgumentParser()
parser.add_argument("-f","--file",required=True)
args = parser.parse_args()

x=list()
y=list()

# read csv file
with open(args.file,"r") as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        x.append(row[0])
        y.append(row[2])

# plot 
plt.plot(x,y,"r-")
plt.show()

