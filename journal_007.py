"""
Created on Wed Jan 19 10:08:47 2022

@author: menglong he 
"""

import numpy as np
from numba import jit
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

X_c = 2
r = 0.189 
mu = 25 
c = 3e10                             # Speed of light in cm / s
gamma = 3.1e-9                       # cm/W
ng = 3.6                             # Refractive index of Silicon at 1.55um. 
vg = c/ng                            # Group velocity
theta = 0.01                         # Assuming a set % of power is coupled to ring from bus 
Lambda_rp = 1550e-7                  # cm 1550 nm 
wr = 2 * np.pi * c / Lambda_rp       # rad/s 1550 nm
alpha_dB = 1 
Radius = 50e-4                       # good choose for critical coupling 
tR = (2 * np.pi * Radius) / vg
Aeff = 1e-9
wr = 2 * np.pi * c / Lambda_rp  
fsr = Lambda_rp**2/(ng*2*np.pi*Radius)*1e7

def threshold (Radius1, Alpha1):
    alpha1 = np.log(10)/10*Alpha1
    a = np.exp(-1*alpha1*2*np.pi*Radius1)
    theta = 1 - a 
    return theta

loss = np.linspace(0.1, 10, 1000)

y1 =  np.linspace(0, 10, 1000)
y2 =  np.linspace(0, 10, 1000)
y3 =  np.linspace(0, 10, 1000)
y4 =  np.linspace(0, 10, 1000)
y5 =  np.linspace(0, 10, 1000)

for i in range(0,len(loss)):
    y1[i] = threshold (20e-4, loss[i])
    y2[i] = threshold (40e-4, loss[i])
    y3[i] = threshold (60e-4, loss[i])
    y4[i] = threshold (80e-4, loss[i])
    y5[i] = threshold (100e-4, loss[i])

plt.plot(loss,y1,label= r"Radius R = 20 $\mu$ m")
plt.plot(loss,y2,label= r"Radius R = 40 $\mu$ m")
plt.plot(loss,y3,label= r"Radius R = 60 $\mu$ m")
plt.plot(loss,y4,label= r"Radius R = 80 $\mu$ m")
plt.plot(loss,y5,label= r"Radius R = 100 $\mu$ m")
plt.legend()
plt.ylabel(r"Coupling coefficient $\theta$")
plt.xlabel(r"Linear loss $\alpha$ (dB/cm)")
plt.xscale("log")
plt.yscale("log")






