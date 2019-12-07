#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse as ae
import seaborn as sns
sns.set_style("darkgrid")
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

parser = ae.ArgumentParser()
parser.add_argument("-d","--data",default="./sampling_log")
args = parser.parse_args()

file = open(args.data,"r")
lines = file.readlines()
file.close()

mu = list()
sgm = list()
loss = list()
var = list()

for i in range(len(lines)):
    a,b,c,d = lines[i].split("\t")
    mu.append(float(a))
    sgm.append(float(b))
    loss.append(float(c))
    var.append(float(d))

fig = plt.figure()
ax = Axes3D(fig)

ax.set_xlabel("mu")
ax.set_ylabel("sigma")
ax.set_zlabel("loss")
#ax.set_zscale("log")
#ax.set_zlim(zmax=1.0e+8)

# search min
loss = np.array(loss)
imin = loss.argmin()
print("loss=",loss.min()," mu=",mu[imin]," sgm=",sgm[imin])

ax.plot(mu,sgm,loss,marker="o",linestyle='None')
plt.show()
