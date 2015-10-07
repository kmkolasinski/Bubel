#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import csv


file = "lattice.dat"
#ax = plt.gca(projection='3d')

pscale=1.0
lscale=10.0
fig, ax = plt. subplots()
ax.set_aspect('equal')
desired=[1,2]
with open(file, 'r') as fin:
    reader=csv.reader(fin)
    result=[[(s) for s in row] for i,row in enumerate(reader) if i in desired]

minCorner = map(float,result[0][0].split())
maxCorner = map(float,result[1][0].split())
xWidth = abs(minCorner[0]-maxCorner[0])
yWidth = abs(minCorner[1]-maxCorner[1])
zWidth = abs(minCorner[2]-maxCorner[2])

ax.set_xlim([minCorner[0]-xWidth*0.1,maxCorner[0]+xWidth*0.1])
ax.set_ylim([minCorner[1]-yWidth*0.1,maxCorner[1]+yWidth*0.1])

data = np.loadtxt(file,skiprows=4)
no_lines = np.size(data[:,0])
#no_lines = 200
for i in range(no_lines):
    x = [data[i,0],data[i,3]]
    y = [data[i,1],data[i,4]]
    #z = [data[i,2],data[i,5]]
    ax.plot(x,y,c='k',lw=data[i,6]*lscale,ls='-')              

for i in range(no_lines):  
    if(data[i,6] > 1.0):
        ax.plot(data[i,0],data[i,1],c='r',markersize=data[i,6]*pscale,marker='.',markeredgecolor='k')                  

plt.xlabel("x")
plt.xlabel("y")
plt.savefig("lattice.png")        
#plt.show()        