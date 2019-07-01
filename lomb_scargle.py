#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 11:23:54 2018

@author: jishnu
"""
#%%
import numpy as np
from astropy.stats import LombScargle
import matplotlib.pyplot as plt

#%%
data=np.loadtxt(
        '/home/jishnu/path_to_file/',
        delimiter='\t',
        skiprows=1
        ).astype(np.float)

t=data[:,0]
y=data[:,1]
e=data[:,2]

#%%
frequency,period,periody,power = [],[],[],[]
Period = 0

#%%
def lomb_scargle(t,y):
    global frequency,power,period,periody,Period
    frequency, power = LombScargle(t, y).autopower(maximum_frequency=8)
    period = 1/frequency
    periody = period/365
    if 1==1:
        plt.figure(figsize = (9,6))
        plt.xlabel('Frequency')
        plt.ylabel('Power')
        #plt.xlim(1,6)
        plt.plot(frequency, power)
        plt.show()
        #
    Period = 1/frequency[np.argmax(power)]
    return 1/frequency[np.argmax(power)]

#%%
def phase_fold(t,y,period):
    phases=np.remainder(t,period)/period
    phases=np.concatenate((phases,phases+1))
    y=np.concatenate((y,y))
    plt.figure(figsize = (9,6))
    plt.title('Period : %.6f'%period)
    plt.xlabel('Phase')
    plt.ylabel('Flux')
    plt.grid()
    #plt.ylim(-0.5,0.2)
    plt.gca().invert_yaxis()
    plt.scatter(phases,y,c='k',lw=.1)
    plt.show()
    plt.close()

#%%
def plot_best_period(t,y):
    phase_fold(t,y,lomb_scargle(t,y))

#%%
plot_best_period(t,y)
print('period:',lomb_scargle(t,y))
