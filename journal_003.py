# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:01:40 2022

@author: mengl
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:05:03 2022

@author: mengl
"""

# Importing packages
import matplotlib.pyplot as plt
import numpy as np
import math

r = 0.2
mu = 25
Aeff = 1e-9
c = 3e10 # cm/s
ng = 3.6 
vg = c/ng 
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
def n(radius1,loss,theta1):
    tR = (2 * np.pi * radius1) / vg # s 
    tau_ph1 = photonlifetime(radius1,theta1,loss) 
    print(tau_ph1*1e12)
    normal_f1 = (2*gamma*vg*tau_ph1)**(-0.5)
    normal_in = ((tR**2)/(8*theta*gamma*vg*(tau_ph1**3)))**(0.5)
    normal_tau = 2*tau_ph1 
    return (normal_f1**2)*Aeff

def si(xc):
    xc = xc * 5  
    a = -1*xc/mu
    b = (1 + (r**2))**0.5 - (2*r)
    c = -1 
    delta = (b**2) - (4*a*c)
    ap1 = ((-1*b) + np.sqrt(delta))/(2*a)
    ap2 = ((-1*b) - np.sqrt(delta))/(2*a)
    return ap1,ap2


x= np.linspace(1e-2,1e0,2000)
y12,y11 = si(x)
y11 = y12 * n(20e-4,0.7,0.1)
y12 = y12 * n(50e-4,0.7,0.1)
y13 = y12 * n(200e-4,0.7,0.1)

plt.plot(x, y11,"b",label=(r"Power range with R = 20 $\mu$m"))
plt.plot(x, y12,"r",label=(r"Power range with R = 50 $\mu$m"))
plt.plot(x, y13,"g",label=(r"Power range with R = 100 $\mu$m"))
# plt.plot(x, y12,"r",label=("Lower border of the parametric gain region"))
plt.xlabel(r"Carrier lifetime / photon lifetime ($\tau{c}$/$\tau{ph}$)")
plt.ylabel("Circulating power $P_{c}$ (W)")
plt.xscale("log")
plt.legend()
plt.show()