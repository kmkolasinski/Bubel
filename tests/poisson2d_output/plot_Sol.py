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

filename = "Potential.dat"

data     = np.loadtxt(filename)
plt.clf()
x     = data[:,0]
y     = data[:,1]

no_x_samples = 100
no_y_samples = 100

xmin , ymin = min(x) , min(y)
xmax , ymax = max(x) , max(y)

xi = np.arange(xmin, xmax , (xmax-xmin)/no_x_samples)
yi = np.arange(ymin, ymax , (ymax-ymin)/no_y_samples)

zi = griddata(x, y, data[:,2], xi, yi)
ax = plt.subplot(121)
ax.contourf(xi, yi, zi, 80 )
ax.set_aspect('equal')
ax.set_title("Input density $ \\rho(x,y)$")
ax.set_xlabel("x [nm]")
ax.set_ylabel("y [nm]")

zi = griddata(x, y, data[:,3], xi, yi)
ax = plt.subplot(122)
ax.contourf(xi, yi, zi, 80 )
ax.set_aspect('equal')
ax.set_xlabel("x [nm]")
ax.set_title("Solution of $\\nabla^2 \psi = \\rho$ ")

plt.savefig("Potential.png")
#plt.show()
