# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:05:03 2022

@author: mengl
"""

# Importing packages
import matplotlib.pyplot as plt
import numpy as np
X_c = 1
r = 0.189 
mu = 25 

def si(ap):
    return -1 - (X_c*(ap**2)/mu) - (2*r*ap) + ((1+(r**2))**0.5*ap)

def sin(ap):
    return ap - 1 

x= np.linspace(0,30,2000)
y1 = si(x)
y2 = sin(x)
plt.plot(x, y1,"r",label=("Gain of Silicon Ring Cavities"))
plt.plot(x, y2,"b",label=("Gain of Silicon Nitride Ring Cavities"))
plt.xlabel("Normolized circulating power")
plt.ylabel("Gain coefficient G")

plt.legend()
plt.show()