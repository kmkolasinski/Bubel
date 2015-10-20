#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv


file = "bands.dat"

data = np.loadtxt(file)
no_lines = np.size(data[0,:])
x = data[:,0]

ax  =plt.subplot(121)
for i in range(no_lines-1):    
    ax.plot(x,data[:,i+1],c='k',ls='-')     
    
ax.set_xlabel("k [1/unit size]")    
ax.set_ylabel("Energy [some units]")   
ax.set_ylim([-3,3])


file = "T.dat"

data = np.loadtxt(file)
no_lines = np.size(data[0,:])
x = data[:,0]

ax  =plt.subplot(122 , sharey=ax)
ax.plot(data[:,1],x,c='k',ls='-')     
ax.set_ylim([-3,3])
plt.savefig("bands.png")
#plt.show()        