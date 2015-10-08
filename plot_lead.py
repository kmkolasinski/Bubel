#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv
from xml.dom import minidom
from matplotlib.collections import LineCollection

xmldoc = minidom.parse('lead.dat')

file = "lattice.dat"
xml_shape_type  = xmldoc.getElementsByTagName('shape_type')
xml_shape_data  = xmldoc.getElementsByTagName('shape_data')
xml_lead_vector = xmldoc.getElementsByTagName('lead_vector')

print "Shape type:",xml_shape_type[0].childNodes[0].nodeValue

shape_data = xml_shape_data[0].childNodes[0].nodeValue
f_shape_data = shape_data.split()
f_shape_data = [float(i) for i in f_shape_data]
print  f_shape_data


xml_lead_vector = xml_lead_vector[0].childNodes[0].nodeValue
f_lead_vector   = xml_lead_vector.split()
f_lead_vector   = [float(i) for i in f_lead_vector]
print f_lead_vector

clf()
f = plt.figure(1)    
ax  = plt.subplot(111)   
shape_contour_x =  [f_shape_data[0],f_shape_data[2],f_shape_data[2],f_shape_data[0],f_shape_data[0]] 
shape_contour_y =  [f_shape_data[1],f_shape_data[1],f_shape_data[3],f_shape_data[3],f_shape_data[1]] 
ax.plot( shape_contour_x , shape_contour_y )
offset_x = [f_lead_vector[0]]*5
offset_y = [f_lead_vector[1]]*5
ax.plot( map(add, shape_contour_x, offset_x)  , map(add, shape_contour_y, offset_y) )


# rebuild ends using none to separate line segments
xlist = []
ylist = []
lead_datas = xmldoc.getElementsByTagName('lead_data')
for ldatas in lead_datas:
    ldata = ldatas.getElementsByTagName('data')
    print len(ldata)
    for i in range(len(ldata)):        
#    for i in range(6):                
        fdata = ldata[i].childNodes[0].nodeValue.split()
        fdata = fdata[0:7]
        fdata = [float(j) for j in fdata]
        print i,fdata[0:6]
        xlist.extend([fdata[0],fdata[3]])        
        xlist.append(None)
        ylist.extend([fdata[1],fdata[4]])
        ylist.append(None)    

del xlist[-1]
del ylist[-1]
    
ax.plot(xlist,ylist,'b-')



# rebuild ends using none to separate line segments
wlist = []
lines = []
lead_datas = xmldoc.getElementsByTagName('next_cell_lead_data')
for ldatas in lead_datas:
    ldata = ldatas.getElementsByTagName('data')
    print len(ldata)
    for i in range(len(ldata)):        
        fdata = ldata[i].childNodes[0].nodeValue.split()
        fdata = fdata[0:7]
        fdata = [float(j) for j in fdata]
        print i,fdata[6]
        lines.append([ (fdata[0],fdata[1]) , (fdata[3],fdata[4]) ])
        wlist.extend([fdata[6]])


lc = LineCollection(lines, linewidths=wlist,colors='green')
ax.add_collection(lc)


# rebuild ends using none to separate line segments
wlist = []
lines = []
lead_datas = xmldoc.getElementsByTagName('lead_coupling')
for ldatas in lead_datas:
    ldata = ldatas.getElementsByTagName('data')
    print len(ldata)
    for i in range(len(ldata)):        
        fdata = ldata[i].childNodes[0].nodeValue.split()
        fdata = fdata[0:7]
        fdata = [float(j) for j in fdata]
        print i,fdata[6]
        lines.append([ (fdata[0],fdata[1]) , (fdata[3],fdata[4]) ])
        wlist.extend([fdata[6]])


lc = LineCollection(lines, linewidths=wlist,colors='red')
ax.add_collection(lc)

ax.set_aspect('equal') 
ax.margins(0.1)
print(len(xml_lead_data))
