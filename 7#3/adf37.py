#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 20:05:02 2022

@author: cg0527
"""

import os
import time
import numpy as np
import adas as adas
import time
import dist

def integral(x,y):#Calculating the integral value by accumulating the product of the interval length and the function value at the middle of the interval.
    l=len(x)
    I=0
    for i in range(0,l-1):
        I=I+((y[i]+y[i+1])/2.0)*(x[i+1]-x[i])
    return I

def M1(T,x):# Maxwellian diastribution
    me=9.1e-31#mass of electron
    x=np.array(x)
    C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)
    #C=2.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)
    dis=C*x*(1.6e-19)*np.exp(-x/T)
    return np.array(dis)/np.sqrt((2*x*(1.6e-19))/me)


def M(T,x):# Maxwellian diastribution
    me=9.1e-31#mass of electron
    x=np.array(x)
    x1=[]
    #C=4.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)*(1/(2*me))**(1/2)*100.0
    C=2.0*np.pi*(1/(np.pi*T*1.6e-19))**(3/2.0)
    dis=C*(x*(1.6e-19))**0.5*np.exp(-x/T)
    return dis

def nor(mu,s,x):
    x1=(x-mu)/s
    temp=np.exp(-(x1**2.0)/2.0)/np.sqrt(2.0*np.pi)
    return (1/(s*1.6e-19))*temp

def newdis(mu,s,T,r,x):
    return (r*nor(mu,s,x)+(1-r)*M(T,x))*1.6e-19

def transpose(array):#new function for transposing 2d array
    l1=len(array)
    l2=len(array[0])
    NewArray=np.zeros((l2,l1))
    for i in range(0,l1):
        for j in range(0,l2):
            NewArray[j][i]=array[i][j]
    return NewArray


def noNan(array):
    for i in range(len(array)):
        for j in range(len(array[0])):
            if(np.isnan(array[i][j])):
                array[i][j]=0
    return array

def integral(x,y):#Calculating the integral value by accumulating the product of the interval length and the function value at the middle of the interval.
    l=len(x)
    I=0
    for i in range(0,l-1):
        I=I+((y[i]+y[i+1])/2.0)*(x[i+1]-x[i])
    return I

def make_adf37(e,eedf):
    #energy=np.linspace(0.1,1e4,1000)
    h=6.626*10**(-34)
    v_c=3*10**(8)
    Needf=[]
    Ne=[]
    for i in range(len(eedf)):
        if(eedf[i] > 1*10**(-100)):
            Ne.append(e[i])
            Needf.append(eedf[i])
    Ne=np.array(Ne)
    Needf=np.array(Needf)
    nenergy=len(Ne)
    fulldata={'title' : "test Maxwellian", 'icateg' : 2, 'ieunit': 2,
    'nenerg':nenergy, 'nblock': 1, 'nform1': 1, 'nform2': 3, 'param2': 1.0,
    'ea' : Ne, 'fa' : Needf}
    my_file='/home/cg0527/7#3/test_write_adf37.dat'
    adas.write_adf37(file=my_file, fulldata=fulldata)

def rline_read(file):
    f=open(file,'r')
    rate=[]
    with open(file) as f:
        temp = f.read().splitlines()
    for i in range(0,len(temp)):
        if(temp[i][0]=='R'):
            temp1 = [temp[i][2]+temp[i][3],temp[i][6]+temp[i][7],temp[i][17]+temp[i][18]+temp[i][19]+temp[i][20]+temp[i][21]+temp[i][22]+temp[i][23]]
            rate.append(temp1)

    for i in range(0,len(rate)):
        power=0
        base=0
        if(rate[i][2][4]+rate[i][2][5]+rate[i][2][6]=='NaN'):
            print("NaN issue")

        if(rate[i][2][4]+rate[i][2][5]+rate[i][2][6]!='NaN'):
            power=int(rate[i][2][4]+rate[i][2][5]+rate[i][2][6])
            base=int(rate[i][2][0]+rate[i][2][2]+rate[i][2][3])/10
        lev=int(rate[i][0])
        ion=int(rate[i][1])
        rate[i]=[lev,ion,base*10**(power)]
    totalRate=0
    for i in range(len(rate)):
        totalRate=totalRate+rate[i][2]
    return rate,totalRate

#%%
k=4
z=6
TeIndex=10
TeNum=35
r=0
mu=200
dens=1e13
delta_E=198310.7722/8066
s=10
c=0
te=np.linspace(6,20,14)
e=np.linspace(0.1,5e3,300)
ScaledForRate=[]
ForRate=[]
T=te[3]
dis=newdis(mu,s,T,r,e)
dis=dis/integral(e,dis)
print(integral(e,dis))
make_adf37(e,dis)
os.system("sh /home/cg0527/7#3/pass.sh")
#%%
temp=adas.xxdata_04("/home/cg0527/Documents/Urop/test/c0_kappa26.dat")
dict5=adas.xxdata_04("/home/cg0527/7#3/he_rate.dat")
r=dict5["ion"][0][0]
ScaledForRate.append(r)
ForRate.append(r/np.exp(delta_E/T)) 
    
rate,totalRate=rline_read(file="rlines.txt")
parent=[]
lvl=[]
rec=[]
for j in range(len(rate)):
    parent.append(int(rate[j][1]))
    lvl.append(int(rate[j][0]))
    rec.append([rate[j][2]])



#%%
'''
dict5["rec"]=np.array(rec)
dict5["level_rec"]=lvl
dict5["parent_rec"]=parent
dict5["adf04_type"]=4
file="/home/cg0527/7#3/test_write_adf04.dat"
adas.write_adf04(file=file, fulldata=dict5)
temp=adas.xxdata_04("/home/cg0527/Documents/Urop/test/c0_kappa26.dat")

data=adas.xxdata_04(file)
tempte=data['te']
tempte=tempte[0]

file2="he_rate.dat"
adas.run_adas218(adf04=file2,te=9.53e4,dens=dens,pass_dir="/home/cg0527/a211",gcr=1)
'''
     