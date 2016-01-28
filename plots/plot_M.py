#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""


import numpy as np
import pylab

import matplotlib.pyplot as plt
import sys , os
import matplotlib
from matplotlib.mlab import griddata
from scipy.linalg import eig
from scipy.linalg import schur
from scipy import linalg

def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]




def read_cmpx_data(file):
    data = np.loadtxt(file)
    nr = size(data,0)
    nc = size(data,1)
    print "reading mat:",nr,nc
    mat = np.zeros([nr,nc/2],dtype=complex)
    for i in range(nr):
        for j in range(nc/2):
            mat[i,j] = complex(data[i,j*2],data[i,j*2+1])
    return mat

#==============================================================================
# 
# w = read_cmpx_data("fort.110")
# U = read_cmpx_data("fort.111")
# fE =read_cmpx_data("fort.112")
# 
# Linv = inv(diag(w[:,0]))
# Uinv = inv(U)
# 
# E = U.dot(Linv.dot(Uinv))
# print "cond=",cond(U)
# print "Roznica1:",abs(sum(abs(E-fE)))
# 
# w   = read_cmpx_data("fort.310")
# U   = read_cmpx_data("fort.311")
# U   = transpose(U)
# 
# u,s,vt = svd(U)
# Linv = inv(diag(w[:,0]))
# Ilambda = vt.dot(Linv.dot(transpose(conjugate(vt))))
# n = size(s)
# eps = 1.0e-30
# m = size(s[s>eps])
# #m = n
# Ilambda = Ilambda[:m,:m]
# S   = diag(s[s>eps])
# 
# 
# SlS = S.dot(Ilambda.dot(inv(S)))
# 
# L = np.zeros([n,n],dtype=complex)
# L[:m,:m] = SlS
# 
# E2 = u.dot(L.dot(transpose(conjugate(u))))
# 
# print "cond=",cond(L)
# print "Roznica2:",abs(sum(abs(E2-fE)))
# 
# 
# fE2 =read_cmpx_data("fort.312")
# print "Roznica3:",abs(sum(abs(fE2-fE)))
#==============================================================================
#==============================================================================
# a=abs(Ilambda)
# a0 = abs(SlS)
# a1= abs(fE2)
# a2= abs(E2)
#==============================================================================


x    = read_cmpx_data("fort.222")

plt.matshow(abs(x))
w    = read_cmpx_data("fort.223")
plt.matshow(abs(w))
a = abs(w)


#==============================================================================
# U    = read_cmpx_data("fort.510")
# #UL   = read_cmpx_data("fort.313")
# U    = transpose(U)
# #UL   = transpose(UL)
# 
# u,s,vt = svd(U)
# #ul,sl,vlt = svd(UL)
# print cond(u),cond(vt)
# #print cond(ul),cond(vlt)
# 
# #ULinv = inv(UL)
# #E = U.dot(ULinv)
# #print "cond=",cond(U),cond(ULinv)
# 
# #fE   = read_cmpx_data("fort.312")
# #print "Roznica1:",abs(sum(abs(E-fE)))
# 
# #p=fE.dot(U) - U.dot(inv(diag(w[:,0])))
# 
# 
# n = size(U,0)
# Linv = inv(diag(w[:,0]))
# 
# vlv = vt.dot(Linv.dot(transpose(conjugate(vt))))
# print "vlv=",sum(vlv)    
# print "u=",sum(u)
# print "vt=",sum(vt)
# 
# S = diag(s)
# sPs = S.dot(vlv.dot(inv(S)))
# print "sps=",sum(sPs)
# 
#==============================================================================
#==============================================================================
# E = u.dot(sPs.dot(transpose(conjugate(u))))
# print "E=",sum(E)
#==============================================================================
