#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy.linalg import inv
plt.cla()
file = "bands.dat"

data = np.loadtxt(file)
no_lines = np.size(data[0,:])
x = data[:,0]
for i in range(no_lines-1):
    plt.plot(x,data[:,i+1],c='k',ls='-')     
    
x = np.array(x)
#plt.plot(x,0.2+x*0)    
#plt.xlim([-0.5,0.5])
#plt.ylim([0,0.2])
plt.xlabel("k [1/unit size]")    
plt.ylabel("Energy [some units]")   
plt.savefig("bands.pdf")
#plt.show()        


#data = np.loadtxt("fort.2")
#n = np.size(data[:,0])
#
#cdata = np.zeros([n,n],dtype=np.complex64)
#cphi = np.array([0]*n,dtype=np.complex64)
#for i in range(n):
#    cphi[i] = np.complex( data[i,0] , data[i,1])
#    for j in range(n):
#        cdata[i,j] = np.complex( data[i,2+2*j] , (data[i,2+2*j+1]))
#    
#    
#ainv = inv(np.matrix(cdata))    
#
#c =ainv.dot(cphi)
#for i in range(n):
#    c[0,i] = np.abs(c[0,i])**2
#print abs(c)