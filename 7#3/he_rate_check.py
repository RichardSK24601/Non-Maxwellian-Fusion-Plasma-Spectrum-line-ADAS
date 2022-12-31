#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:01:21 2022

@author: cg0527
"""
import numpy as np
import matplotlib.pyplot as plt
import adas
from mpmath import *#Package for high precision calculation
mp.dps = 25; mp.pretty = True#setting up mpmath package
import scipy as sp


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


def conv(x,y1,y2):# convolution function
    y1=np.array(y1)
    y2=np.array(y2)
    temp=y1*y2
    return integral(x, temp)
        
def rate(e,c,dis):# caculate the excitation rate
    me=9.1e-31
    e=np.array(e)
    c=np.array(c)*(2*e*1.6e-19/me)**0.5*100
    e=e*1.6e-19
    return conv(e,c,dis)

def integral(x,y):#Calculating the integral value by accumulating the product of the interval length and the function value at the middle of the interval.
    l=len(x)
    I=0
    for i in range(0,l-1):
        I=I+((y[i]+y[i+1])/2.0)*(x[i+1]-x[i])
    return I
def smooth(e1,dis,e2):
    Ndis=[]
    for i in range(0,len(e2)):
        if(e2[i]<e1[0]):
            Ndis.append(0)
        if(e2[i]>=e1[len(e1)-1]):
            Ndis.append(dis[(len(e1)-1)])
        for j in range(0,len(e1)-1):
            if(e2[i]>=e1[j] and e2[i]<e1[j+1]):
                temp=dis[j]+((dis[j+1]-dis[j])/(e1[j+1]-e1[j]))*(e2[i]-e1[j])
                Ndis.append(temp)
    if(len(Ndis)!=len(e2)):
        print("wrong length")
    return np.array(Ndis)
file="/home/adas/adas/adf04/helike/helike_idp04he0_t1.dat"
#file="/home/adas/adas/adf04/adas#2/mom97_ls#he0.dat"
file2="/home/adas/adas/adf11/scd96/scd96_he.dat"
dict1=adas.xxdata_04(file)
const1=13.6
const2=8.7972e-17
delta_E=198310.7722/8066
g_i=1
dens=1e6
X=dict1["te"]
omega=dict1["ion"][0]
X=np.array(X)
omega=np.array(omega)
sigma=(omega/g_i)*const1*const2/(delta_E*X)
e=X*delta_E
Ne=np.linspace(0.1,1e5,50000)
Nsigma=smooth(e,sigma,Ne)
te=np.linspace(6,20,14)
adf07=[]
adf11=[]
pyconv=[]
ScaledPyConv=[]

for i in range(0,len(te)):
    T=te[i]
    dis=M(T,Ne)
    r=rate(Ne,Nsigma,dis)                                                                                                                                                                                                                                                                                                                                         
    dict2=adas.read_adf11(file=file2, adf11type='scd', is1=1, te=T, dens=dens)
    dict3=adas.xxdata_04("/home/adas/adas/adf04/adas#2/mom97_ls#he0.dat")
    dict4=adas.read_adf07("dere07#he.dat",te=T,unit_te='ev')
    adf07.append(dict4[0][0])
    adf11.append(dict2[0])
    pyconv.append(r)
    ScaledPyConv.append(r*np.exp(delta_E/T))
    print('Ionization rate:',r)
    print('Scaled Ionization rate:',r*np.exp(delta_E/T))
    print("ion rate from adf11:",dict2[0])
    print("ion rate from adf07:",dict4[0][0])
    #print("Normalization check:",integral(Ne,M(T,Ne))*1.6e-19)
'''
plt.figure()
plt.plot(te,adf07,label='adf07')
plt.plot(te,adf11,label='adf11')
plt.plot(te,pyconv,label='pyRate')
plt.plot(te,ScaledPyConv,label='ScaledPyRate')
plt.plot(te,ScaledForRate,label='ScaledForRate')
plt.plot(te,ForRate,label='ForRate')
plt.yscale("log")
plt.ylabel("ionization rate/ cm^3s^-1")
plt.xlabel("Temperature/eV")
plt.legend()
'''