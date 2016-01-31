#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys



nTest = sys.argv[1]
file = "Test"+nTest+".dat"
print "Plotting test:",file


plt.clf()
data = np.loadtxt(file)
no_lines = np.size(data[0,:])
x = data[:,0]
T_exact   = data[:,1]
eps = 10e-15
dT_auto   = abs(T_exact-data[:,2]+eps) 
dT_ggev   = abs(T_exact-data[:,3]+eps) 
dT_schur  = abs(T_exact-data[:,4]+eps) 
#plt.plot(x,T_exact,c='k',ls='-')  
#plt.plot(x,data[:,3],c='k',ls='-')  
ax = plt.subplot(211)
plt.yscale("log")   
plt.plot(x,dT_auto,c='r',ls='-',label="Auto  method")        
plt.plot(x,dT_ggev,c='k',ls='-',label="GGEV  method")        
plt.plot(x,dT_schur,c='b',ls='-',label="Schur method")   


ax2 = ax.twinx()
ax2.plot(x, T_exact, 'r.')
ax2.set_ylabel('sin', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

ax.legend(loc="lower right")     
ax.set_ylabel("Error")   
 
plt.ylabel("Num. modes")   

ax = plt.subplot(212)
   
plt.plot(x,data[:,5],c='r',ls='-',label="Auto  method")        
plt.plot(x,data[:,6],c='k',ls='-',label="GGEV  method")        
plt.plot(x,data[:,7],c='b',ls='-',label="Schur method")   
ax.set_xlabel("Ef [meV]")   
ax.set_ylabel("CPU time [s]") 
#ax.legend() 
plt.savefig(file+".png")
#plt.show()        
