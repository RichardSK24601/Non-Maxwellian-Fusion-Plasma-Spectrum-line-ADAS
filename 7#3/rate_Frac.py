#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 11:16:38 2022

@author: cg0527
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
from scipy.linalg import null_space
from mpmath import *#Package for high precision calculation
from scipy.special import gamma
import Ratio
import adas
mp.dps = 25; mp.pretty = True#setting up mpmath package

def Kappa(k,E_k,x):
    me=9.1e-31
    A=(2/(E_k*1.6e-19*np.sqrt(np.pi)))*np.sqrt(x/E_k)
    B=k**(-3/2)*(gamma(k+1))/gamma(k-0.5)
    C=(1+(x/(k*E_k)))**(-k-1)
    D=np.sqrt(2*x*1.6e-19/me)
    return A*B*C#*D*100
def M1(T,x):# Maxwellian diastribution
    me=9.1e-31#mass of electron
    x=np.array(x)
    x1=[]
    C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)*100.0
    #C=2.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)
    for i in range(len(x)):
        dis=C*x[i]*(1.6e-19)*exp(-x[i]/T)
        x1.append(dis)
    return np.array(x1)

def M(T,x):# Maxwellian diastribution
    me=9.1e-31#mass of electron
    x=np.array(x)
    x1=[]
    #C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)*100.0
    C=2.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)
    for i in range(len(x)):
        dis=C*(x[i]*(1.6e-19))**0.5*exp(-x[i]/T)
        x1.append(dis)
    return np.array(x1)
def nor(mu,s,x):
    x1=(x-mu)/s
    temp=np.exp(-(x1**2.0)/2.0)/np.sqrt(2.0*np.pi)
    return (1/(s*1.6e-19))*temp

def newdis(mu,s,T,r,x):
    return r*nor(mu,s,x)+(1-r)*M(T,x)

def integral(x,y):#Calculating the integral value by accumulating the product of the interval length and the function value at the middle of the interval.
    l=len(x)
    I=0
    for i in range(0,l-1):
        I=I+((y[i]+y[i+1])/2.0)*(x[i+1]-x[i])
    return I

def conv(x,y1,y2):# convolution function
    y1=np.array(y1)
    y2=np.array(y2)
    temp=y1*y2
    return integral(x, temp)
        
def rate(e,c,dis):# caculate the excitation rate
    me=9.1e-31
    e=np.array(e)
    D=np.sqrt(2*e*1.6e-19/me)
    c=np.array(c)
    c=c*D*100
    e=e*1.6e-19
    return conv(e,c,dis)

def smooth(e1,dis,e2):
    Ndis=[]
    for i in range(0,len(e2)):
        if(e2[i]<e1[0]):
            Ndis.append(e1[0])
        if(e2[i]>=e1[len(e1)-1]):
            Ndis.append(dis[(len(e1)-1)])
        for j in range(0,len(e1)-1):
            if(e2[i]>=e1[j] and e2[i]<e1[j+1]):
                temp=dis[j]+((dis[j+1]-dis[j])/(e1[j+1]-e1[j]))*(e2[i]-e1[j])
                Ndis.append(temp)
    if(len(Ndis)!=len(e2)):
        print("wrong length")
    return np.array(Ndis)
def transpose(array):#new function for transposing 2d array
    l1=len(array)
    l2=len(array[0])
    NewArray=np.zeros((l2,l1))
    for i in range(0,l1):
        for j in range(0,l2):
            NewArray[j][i]=array[i][j]
    return NewArray
def frac(fileacd,filec,z,dens,dis,disEgrid,te):
    res = np.load(filec, allow_pickle=True)
    E = res['energy']
    C = res['sigma']
    acd=[]
    for i in range(0,z):
        temp=adas.read_adf11(file=fileacd, adf11type='acd', is1=i+1, te=te, dens=dens, all=True)
        acd.append(temp[0][0])
    NewE=transpose(E)
    NewC=transpose(C)
    l3=len(NewE)
    Rate=[]
    temp=[]
    for i in range(0,l3):
        temp1=rate(NewE[i],NewC[i],smooth(disEgrid,dis,NewE[i])) 
        #print(temp1)
        temp.append(temp1)
    Rate=np.array(temp)    
    return Ratio.ratio(acd,Rate)

filec='c_cross_section.npz'
fileacd='/home/adas/adas/adf11/acd96/acd96_c.dat'
dens=1e13
TeNum=100
z=6
k=10
res = np.load(filec, allow_pickle=True)
te=np.geomspace(1*10**(-0.5), 1.0e4, TeNum)
fraction=[]
r=[]
E = res['energy']
Ne=np.linspace(0.1,2e4,10000)
print(integral(Ne,Kappa(k,te[13],Ne))*1.6e-19)


for i in range(0,len(te)):
    T=te[i]
    E = res['energy']
    C = res['sigma']
    NewE=transpose(E)
    NewC=transpose(C)
    dis=M(T,NewE[0])
    I=integral(NewE[0],dis)
    f=frac(fileacd,filec,z,dens,dis,transpose(E)[0],te[i])
    for j in range(0,len(f)):
        if(f[j]<1e-4):
            f[j]=0
    fraction.append(f)

plt.plot(te,fraction)
plt.xscale("log")
plt.yscale("log")
