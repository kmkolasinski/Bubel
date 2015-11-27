#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv


file = "T.dat"
plt.clf()
data = np.loadtxt(file)
no_lines = np.size(data[0,:])
x = data[:,0]
plt.plot(x,data[:,1],c='k',ls='-')     
  
plt.xlabel("Ef [energy units]")    
plt.ylabel("Transmission")   
plt.savefig("T.png")
#plt.show()        

