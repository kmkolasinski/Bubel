#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv

fig = plt.figure(num=1, figsize=(8, 3), dpi=256, facecolor='w', edgecolor='k')

plt.clf()
def plot_band(filename,nr):
    data = np.loadtxt(filename)
    no_lines = np.size(data[0,:])
    x = data[:,0]
    
    ax  =plt.subplot(nr)
    for i in range(no_lines-1):    
        ax.plot(x,data[:,i+1],c='k',ls='-')     
        
    ax.set_xlabel("k [1/unit size]")    
    #ax.set_ylabel("Energy [some units]")   
    ax.set_ylim([-3,3])
    return ax

plot_band("bands_e.dat",131)
plot_band("bands_h.dat",132)
plot_band("bands_s.dat",133)

plt.subplots_adjust(bottom=0.2)
plt.savefig("bands.png")
#plt.show()        