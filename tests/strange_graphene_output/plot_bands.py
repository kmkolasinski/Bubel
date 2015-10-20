#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt


def plot_bands(file):
    data = np.loadtxt(file+".dat")
    no_lines = np.size(data[0,:])
    x = data[:,0]
    for i in range(no_lines-1):
        plt.plot(x,data[:,i+1],c='k',ls='-')     
        
    plt.xlabel("k [1/unit size]")    
    plt.ylabel("Energy [some units]")   
    plt.savefig(file+".png")


plt.cla();
file = "bands"
plot_bands(file)

#plt.show()        