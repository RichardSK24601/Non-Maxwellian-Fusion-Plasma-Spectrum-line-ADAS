#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 19:46:00 2022

@author: cg0527
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy.linalg import null_space
from mpmath import *#Package for high precision calculation
from scipy.special import gamma
mp.dps = 25; mp.pretty = True#setting up mpmath package
import adf37 as a37

def integral(x,y):#Calculating the integral value by accumulating the product of the interval length and the function value at the middle of the interval.
    l=len(x)
    I=0
    for i in range(0,l-1):
        I=I+((y[i]+y[i+1])/2.0)*(x[i+1]-x[i])
    return I

def M(T,x):# Maxwellian diastribution
    me=9.1e-31#mass of electron
    x=np.array(x)
    x1=[]
    #C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)*100.0
    C=2.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)
    dis=C*(x*(1.6e-19))**0.5*np.exp(-x/T)
    return dis

def pert(delta,T,f0,e):
    me=9.1e-31
    C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)*100.0*(1.6e-19)
    me=9.1e-31
    pth=(me*T)**(0.5)
    deltaP=delta*pth
    p0=pth-deltaP
    p1=pth+deltaP
    E0=0.5*p0**2/me
    E1=0.5*p1**2/me
    pert0=[]
    for i in range(0,len(e)):
        if (e[i]<E0 or e[i]>E1):
            pert0.append(0)
        if(e[i]>E0 and e[i]<E1):
            p=(2*me*e[i])**(0.5)
            temp=np.sin(np.pi*(p-p0)/deltaP) 
            pert0.append(temp)
    pert0=f0*np.array(pert0)*C
    temp= M(T,e)-pert0
    temp=temp/integral(e,temp)
    #for i in range(0,len(temp)):
        #temp[i]=np.format_float_scientific(temp[i])
    return temp

def transpose(array):#new function for transposing 2d array
    l1=len(array)
    l2=len(array[0])
    NewArray=np.zeros((l2,l1))
    for i in range(0,l1):
        for j in range(0,l2):
            NewArray[j][i]=array[i][j]
    return NewArray