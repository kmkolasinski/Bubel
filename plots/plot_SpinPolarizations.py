#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.colors as colors
import matplotlib.cm as cmx

file     = "polarizations.dat"

f = open(file)
lines = f.readlines()
plt.cla()
cmap = plt.cm.jet
scalarMap = cmx.ScalarMappable(cmap=cmap)

for i in range(np.size(lines)):
    data = [float(x) for x in lines[i].split()]
    no_modes = int(data[0])
    Ef    = data[1]           # Fermi energy      
    kvecs = data[2:no_modes+2]# take wave vectors
    pXYZ  = data[no_modes+2:] # take polarizations
    p =   np.reshape(pXYZ,[3,no_modes])
    print "m=",no_modes,"Ef=",Ef
    #print p
    for m in range(no_modes):
         p[:,m] /=  sqrt(sum(p[:,m]**2))   
         
    colors = cmx.jet((p[2,:]+1)/2)
    
    plt.scatter(kvecs,[Ef]*no_modes,c=colors,zorder=2)
    
    for m in range(no_modes):
        #print "m=",m," p=", sqrt(sum(p[:,m]**2))
        c = cmx.jet((p[2,m]+1)/2)
        p[0,m] /= 25.0
        p[:,m] /= 2.0        
        plt.arrow(kvecs[m], Ef, p[0,m], p[1,m], head_width=0.001, head_length=0.01,color=c)
       
       
file = "bands.dat"

data = np.loadtxt(file)
no_lines = np.size(data[0,:])
E0= 27211.384523
x = data[:,0]
for i in range(no_lines-1):
    plt.plot(x,data[:,i+1]*E0,c='k',ls='-',zorder=0)            
   
plt.xlim([-0.4,0.4])
plt.ylim([-0.1,6])
plt.xlabel("k [1/unit size]")    
plt.ylabel("Energy [meV]")   