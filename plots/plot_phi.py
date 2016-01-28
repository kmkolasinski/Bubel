#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: K
"""
import numpy as np
import pylab
import matplotlib.pyplot as plt
import sys , os
import matplotlib
from matplotlib.mlab import griddata

filename = "phi.txt"

data     = np.loadtxt(filename)
plt.clf()
x     = data[:,0]
y     = data[:,1]

dens  = data[:,2]
no_x_samples = 100
no_y_samples = 100

xmin , ymin = min(x) , min(y)
xmax , ymax = max(x) , max(y)

xi = np.arange(xmin, xmax , (xmax-xmin)/no_x_samples)
yi = np.arange(ymin, ymax , (ymax-ymin)/no_y_samples)
zi = griddata(x, y, dens, xi, yi)
ax = plt.subplot(111)
ax.contourf(xi, yi, zi, 100 )
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig("phi.png")
#plt.show()