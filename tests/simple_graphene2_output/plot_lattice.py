#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.collections import LineCollection

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

ax.scatter(minCorner[0],minCorner[1],s=0)
ax.scatter(maxCorner[0],maxCorner[1],s=0)

ax.margins(0.1)
data = np.loadtxt(file,skiprows=4)
no_lines = np.size(data[:,0])
            


        
wlist = []
lines = []
for i in range(no_lines):
        lines.append([ (data[i,0],data[i,1]) , (data[i,3],data[i,4]) ])
        wlist.extend([data[i,6]*lscale])        
        
lc = LineCollection(lines, linewidths=wlist,colors='black',lw=1.0)
ax.add_collection(lc)        


wlist  = []
points = []
for i in range(no_lines):
        if(data[i,6] > 1.0):
            points.append([data[i,0],data[i,1] ])
            wlist.extend([data[i,6]*pscale])        
            
points = np.array(points)
wlist  = np.array(wlist)
if(np.size(points) > 0):        
    ax.scatter(points[:,0],points[:,1], cmap='PuBu', c=wlist , s=50 , edgecolors='k' , zorder=2 )     

plt.savefig("lattice.pdf")