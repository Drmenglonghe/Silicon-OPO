# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 13:03:34 2022

@author: mengl
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
######################################################################
r = 0.2
mu = 25
Aeff = 1e-9
c = 3e10 # cm/s
ng = 3.6 
vg = c/ng 
X_th = 0
Tau_th = 0.04 * 1e-6  
Aeff = 0.1e-8 
neff = Aeff
theta = 1e-1# critical coupling 
radius = 20e-4 # um in cm
alpha_dB = 0.7 # dB/cm
gamma = 3.1e-9

######################################################################
def Linear_response(r1,alpha_dB1):
    alpha1 = np.log(10)/10*alpha_dB1
    return alpha1

def photonlifetime(r2,theta2,alpha_dB1):
    r2 = r2
    tR   = (2 * np.pi * r2) / vg
    L = 2*np.pi*r2
    Loss = Linear_response(r2,alpha_dB1)*L + theta2
    return  tR/Loss

######################################################################
tR = (2 * np.pi * radius) / vg # s 
tau_ph = photonlifetime(radius,theta,alpha_dB) # s 
normal_f = (2*gamma*vg*tau_ph)**(-0.5)
normal_in = ((tR**2)/(8*theta*gamma*vg*(tau_ph**3)))**(0.5)
normal_tau = 2*tau_ph  


def resonance (X_c1,S):   
    S = S/((normal_f**2)*Aeff)
    a = r 
    G = -1 - (X_c1*(S**2)/mu) - (2*r*S) + ((1+(r**2))**0.5*S)
    if G > 0:
        return 0
    else:
        return 1 
##########################################################################################################################

y = np.linspace(1e-1, 10, 1000)
x = np.linspace(5e-2, 50, 1000)
Z0 = np.zeros((len(y),(len(x))))

for i in range(0,len(y)):
    for j in range(0,len(x)):
        Z0[i][j] = resonance(x[j],y[i])

fig, ax = plt.subplots()
x = x/5
im = ax.contourf(x, y, Z0, 36,cmap = "RdGy",linewidths= 0.6)
plt.xscale("log")
plt.yscale("log")
plt.ylabel(r"Circulating pump power $P_{c}$ (W)")
plt.xlabel(r"Free carrier/photon lifetime ($\tau{c}$/$\tau{ph}$)")
plt.title(r" Power threshold for the optimal parametric gain")
#divider = make_axes_locatable(ax)
#cax = divider.new_vertical(size = '4%', pad = 0.5)
#fig.add_axes(cax)
#cbar = fig.colorbar(im, cax = cax, orientation = 'horizontal')
plt.savefig("dopo001.png",dpi=1000,bbox_inches='tight')
plt.show()