#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:07:34 2022

@author: cg0527
"""

import numpy as np

def smooth(e,dis,x):
    if (x<e[0]):
        return dis[0]
    if (x>=e[len(e)-1]):
        return dis([len(e)-1])
    for i in range(0,len(dis)-2):
        if(x>=e[i] and x<e[i+1]):
            temp=dis[i]+((dis[i+1]-dis[i])/(e[i+1]-e[i]))*(x-e[i])
            return temp
        

def getGrid(dis1,e1,dis2,e2,step,maxgirdnum):
    Newdis=[]
    NewE=[]
    if(dis1[0]<=dis2[0]):
        x0=dis1[0]
    else:
        x0=dis2[0]
    if(dis1[len(dis1)-1]<=dis2[len(dis2)-1]):
        xfin=dis2[0]
    else:
        xfin=dis1[0]
    tempE=np.linspace(x0,xfin,maxgirdnum)#much finer than origianl grid