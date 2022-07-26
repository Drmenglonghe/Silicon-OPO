# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:29:55 2022

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
radius = 100e-4 # um in cm
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


def resonance (S,y,theta1):   
    X_c1 = 5*y
    tau_ph = photonlifetime(radius,theta1,alpha_dB) # s 
    normal_f = (2*gamma*vg*tau_ph)**(-0.5)
    normal_in = ((tR**2)/(8*theta1*gamma*vg*(tau_ph**3)))**(0.5)
    normal_tau = 2*tau_ph  
    S = S/((normal_in**2)*Aeff)
    a = r 
    b = X_c1/mu
    p6 = b**2
    p5 = 2*a*b
    p4 = (a**2) + (2*b)
    p3 = (2*a)
    p2 = (1)
    p1 = -(S)
    p = [p6,p5,p4,p3,p2,p1]
    Solution0 = np.roots(p)
    Solution = np.array([ num for num in Solution0 if np.angle(num) == 0 ])
    y = 10*np.log10(max(Solution)*(normal_f**2)*Aeff*1e3)
    return y
##########################################################################################################################

x = np.linspace(1e-3, 1, 100)
y = np.linspace(1e-4, 1e0, 100)
Z0 = np.zeros((len(y),(len(x))))

for i in range(0,len(y)):
    for j in range(0,len(x)):
        Z0[i][j] = resonance(x[j],0.1,y[i])

fig, ax = plt.subplots()
x = 1000 * x
y = y
im = ax.contourf(x, y, Z0, 24,cmap = "gnuplot",linewidths= 0.6)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Input pump power $P_{in}$ (mW)")
plt.ylabel(r"Coupling coefficient $\theta$")
plt.title(r" Circulating pump power at resonance $P_{c}$ (W)")
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size = '4%', pad = 0.5)
fig.add_axes(cax)
cbar = fig.colorbar(im, cax = cax, orientation = 'horizontal')
plt.savefig("dopo001.png",dpi=1000,bbox_inches='tight')
plt.show()