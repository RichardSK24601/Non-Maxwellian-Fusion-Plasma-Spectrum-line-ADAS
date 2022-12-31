#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 11:24:46 2022

@author: cg0527
"""

import numpy as np
import adas
import matplotlib.pyplot as plt

file='/home/adas/adas/adf08/rrc96#h/rrc96#h_he1ls.dat'
adas.run_adas211(file=file, rlines='rlines.txt',adf37="test_write_adf37.dat")
te=[1.16E+04,  2.32E+04,  5.80E+04,  1.16E+05,  2.32E+05,  5.80E+05,  1.16E+06,  2.32E+06,  5.80E+06,  1.16E+07,  2.32E+07,  5.80E+07,  1.16E+08,2.32E+08]
a0=[1.66e-13, 1.18e-13, 7.49e-14, 5.21e-14, 3.44e-14, 1.70e-14, 8.72e-15, 4.07e-15, 1.33e-15, 5.47e-16, 2.16e-16, 5.89e-17, 2.18e-17, 8.03e-18]
a3=[2.04e-13, 1.45e-13, 9.17e-14, 6.42e-14, 4.37e-14, 2.39e-14, 1.37e-14, 7.07e-15, 2.58e-15, 1.12e-15, 4.58e-16, 1.33e-16, 5.08e-17, 1.90e-17]
a20=[1.70e-13, 1.21e-13, 7.64e-14, 5.32e-14, 3.52e-14, 1.77e-14, 9.19e-15, 4.34e-15, 1.44e-15, 5.90e-16, 2.32e-16, 6.51e-17, 2.43e-17, 8.94e-18]
d1=[1.66e-13, 1.18e-13, 7.49e-14, 5.21e-14, 3.44e-14, 1.70e-14, 8.74e-15, 4.09e-15, 1.34e-15, 5.45e-16, 2.14e-16, 5.96e-17, 2.22e-17, 8.15e-18]
d100=[1.41e-13, 9.85e-14, 6.15e-14, 4.37e-14, 2.84e-14, 1.17e-14, 5.08e-15, 2.09e-15, 6.02e-16, 2.28e-16, 8.48e-17, 2.25e-17, 8.10e-18, 2.91e-18]
plt.figure()
plt.plot(te,a0,label='Maxwellian')
plt.plot(te,a3,label='Kappa=3')
plt.plot(te,a20,label='Kappa=20')
plt.plot(te,d1,label='Druyvestyn, x=1')
plt.plot(te,d100,label='Druyvestyn, x=100')
plt.ylabel("recombination rate/ cm^3s^-1")
plt.xlabel("Temperature/eV")
plt.yscale('log')
plt.legend()