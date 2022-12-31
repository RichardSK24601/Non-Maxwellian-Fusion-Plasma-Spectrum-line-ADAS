#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:29:47 2022

@author: cg0527
"""

import os
import subprocess
import pexpect as pex
import time
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy.linalg import null_space
import adas as adas
import numpy as np

####This is where we made all distribution####
###the energy range is expected to be 0-20000eV#
def M(T,x):# Maxwellian diastribution
    me=9.1e-31#mass of electron
    x=np.array(x)
    x1=[]
    C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)
    #C=2.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)
    dis=C*x*(1.6e-19)*np.exp(-x/T)
    return np.array(dis)/np.sqrt((2*x*(1.6e-19))/me)

def nor(mu,s,x):
    x1=(x-mu)/s
    temp=np.exp(-(x1**2.0)/2.0)/np.sqrt(2.0*np.pi)
    return (1/(s*1.6e-19))*temp
def newdis(mu,s,T,r,x):#Max-Nor distribution
    return (r*nor(mu,s,x)+(1-r)*M(T,x))*1.6e-19
def integral(x,y):#Calculating the integral value by accumulating the product of the interval length and the function value at the middle of the interval.
    l=len(x)
    I=0
    for i in range(0,l-1):
        I=I+((y[i]+y[i+1])/2.0)*(x[i+1]-x[i])
    return I

def normal(x,y):
    x=np.array(x)
    y=np.array(y)
    I=integral(x,y)
    return x,y/I

def smooth(E1,f1,E2):
    f2=np.zeros(len(E1))
    for i in range(0,len(E2)):
        for j in range(0,len(E1)-1):
            if(E1[j]<= E2[i] and E2[i]<E1[j+1]):
                f2[i]=f1[j]+(f1[j+1]-f1[j])/(E1[j+1]-E1[j])
    return f2

def transpose(array):#new function for transposing 2d array
    l1=len(array)
    l2=len(array[0])
    NewArray=np.zeros((l2,l1))
    for i in range(0,l1):
        for j in range(0,l2):
            NewArray[j][i]=array[i][j]
    return NewArray


k=4
z=6
TeIndex=10
TeNum=35
r=0.3
mu=200
dens=1e13
s=10
c=0
te=np.geomspace(1*10**(-0.5), 1.0e4, 35)
e=np.linspace(0.1,1e4,1000)
f=newdis(mu,s,te[TeIndex],r,e)
Needf=[]
Needf=[]
Ne=[]
for i in range(len(e)):
    if(f[i] > 1*10**(-30)):
        Ne.append(e[i])
        Needf.append(f[i])
Ne=np.array(Ne)
Needf=np.array(Needf)
nenergy=len(Ne)
fulldata={'title' : "test Maxwellian", 'icateg' : 2, 'ieunit': 2,
'nenerg':nenergy, 'nblock': 1, 'nform1': 1, 'nform2': 3, 'param2': 1.0,
'ea' : Ne, 'fa' : Needf}
my_file='/home/cg0527/7#3/test/test_write_adf37.dat'
adas.write_adf37(file=my_file, fulldata=fulldata)
'''
##random distribution#
filec='be_cross_section.npz'
res = np.load(filec, allow_pickle=True)
E = res['energy']
C = res['sigma']
d=len(E)
E=transpose(E)
Etemp=np.geomspace(10,21000,d)
f=np.random.rand(d)
Etemp,f=normal(Etemp,f)
ftemp=smooth(Etemp,f,E[0])
plt.plot(Etemp,f)
plt.plot(E[0],ftemp)
plt.xscale("log")
plt.yscale("log")
'''