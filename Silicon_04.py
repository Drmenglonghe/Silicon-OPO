# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 13:05:21 2022

@author: menglong He, Germany
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

eta = 1
Aeff = 1e-9
radius = 50e-4
theta = 0.01
alpha_dB = 2
X_c = 0.05
c = 3e10                             
ng = 3.6 
gamma = 3.1e-9
vg = c/ng                         
mua = 25 
r = 0.189
X_th = 0.34
Tau_th = 0.04 * 1e-6  
Lambda_rp = 1550e-7
wr = 2 * np.pi * c / Lambda_rp       
delta_T = 1.8e-4*(wr/c) 


def Linear_response(r1,alpha_dB1):
    alpha1 = 10**(alpha_dB1/10)
    return alpha1 - 1

def photonlifetime(r2,theta2,alpha_dB1):
    r2 = r2
    tR   = (2 * np.pi * r2) / vg
    L = 2*np.pi*r2
    Loss = Linear_response(r2,alpha_dB1)*L + theta2
    return  tR/Loss

tau_ph1 = photonlifetime(radius, theta, alpha_dB)

def silicon_ring_cavity (delta, s, k):
    tR = (2 * np.pi * radius) / vg
    temp = X_th*Tau_th/(2*tau_ph1)
    a = r 
    b = X_c/mua
    c = (X_c -(r*temp))
    d = 2*X_c/mua*temp
    p7 = (d**2)
    p6 = 2*c*d
    p5 = (b**2) + (c**2) + (2*d)
    p4 = (2*a*b) - (2*delta*d) - (2*c)
    p3 = (a**2) + (2*b) + 1 + (2*delta*c)
    p2 = (2*a) - (2*delta)
    p1 = 1 + (delta**2)
    p0 = -(s**2)
    p = [p7,p6,p5,p4,p3,p2,p1,p0]
    
    # a = r
    # b = X_c/mua 
    # c = delta
    # d = -1
    # e = X_c
    # f = ((X_th * Tau_th) /(2*tau_ph1))*r
    # g = ((X_th * Tau_th) /(2*tau_ph1)) * 2* X_c / mua
    # p = [(g**2),(2*f*g)+(2*e*g),(2*e*f)+(2*d*g)+(e**2)+(f**2)+(b**2),(2*d*f)+(2*d*e)+(2*g*c)+(2*a*b),(2*c*f) + (2*e*c)*(d**2)+(2*b)+(a**2),(2*c*d)+(2*a),(c**2) + (1), -s]
    Solution = np.roots(p)
    Solution1 = np.array([ num for num in Solution if np.angle(num) == 0 ])
    flag1 = len(Solution1)
 
    x = max(Solution1)
    phase = delta - (2*x) + (X_c*(x**2))
    D = (1 + (r**2))*(x**2) - ((phase) - (eta*(k**2)))**2
    g = -1 - (2*r*x) - (X_c*(x**2)/mua) + np.sqrt(D)
    g = np.real(g)
    return g

d = 0
Delta_pump = np.linspace(-2,2,100)       # The pump power as control power
Pump1_powers =np.linspace(0,18,100)       # The pump power as control power
z1 = np.empty((len(Pump1_powers),len(Delta_pump)))

for i in range(0,len(Pump1_powers)):
    for j in range(0,len(Delta_pump)):
        z1[i][j] = (silicon_ring_cavity(d,Pump1_powers[i],Delta_pump[j]))

tR = (2 * np.pi * radius) / vg
loss = Linear_response(radius,alpha_dB)*2*np.pi*radius
tau_ph = photonlifetime(radius,theta,alpha_dB)
warning = (8*theta*gamma*vg*(tau_ph**3))
normal_in = ((tR**2)/warning)**(0.5)
Beta_GVD=-0.2*(1e-24/1e2)
normal_t = np.sqrt(abs(Beta_GVD) * vg * tau_ph)

for i in range(0,len(Delta_pump)):
    Delta_pump[i] = (1550e-7 - 2*np.pi*c/(wr - (d / (2*tau_ph)) -(Delta_pump[i]/normal_t)))*1e7
for i in range(0,len(Pump1_powers)):
    Pump1_powers[i] = (Pump1_powers[i]*(normal_in**2)*Aeff*1e3)
    
X, Y = np.meshgrid(Delta_pump,Pump1_powers)
plt.figure(1)
Z1 = z1
 

# plt.contourf(X, Y, Z1,cmap="Reds")
cs1 = plt.contour(X, Y, Z1,alpha=0.6,cmap="Blues",levels = 15)
cs11 = plt.contourf(X, Y, Z1,alpha=0.6,cmap="Blues",levels = 15)
#plt.clabel(cs1, fontsize=10, colors="black")
plt.colorbar()

plt.xlabel(r"Signal pump wavelength difference $\Delta_{\lambda}$ (nm)")
plt.ylabel(r"Pump power $P_{in}$ (mW)")
plt.title(r"Parametric gain in the silicon ring cavities")
plt.savefig("silicon014.png",dpi=1000,bbox_inches='tight')
plt.show()