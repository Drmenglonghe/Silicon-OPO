# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:01:09 2022

@author: mengl
"""


import numpy as np
from numba import jit
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# -*- coding: utf-8 -*-
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

def Linear_response(r1,alpha_dB1):
    alpha1 = np.log(10)/10*alpha_dB1
    return alpha1

def photonlifetime(r2,theta2,alpha_dB1):
    r2 = r2/1e4
    tR   = (2 * np.pi * r2) / vg
    L = 2*np.pi*r2
    Loss = Linear_response(r2,alpha_dB1)*L + theta2
    return  tR/Loss

def Qfactor(r2,alpha,theta2):
    r2 = r2/1e4
    rr2 = (1-theta2)**0.5
    L = 2*np.pi*r2
    a2 = np.exp(-1*Linear_response(r2,alpha)*L)**0.5
    L2 = 2*np.pi*r2
    Q1 = np.pi * ng * L2*np.sqrt(a2*rr2)
    Q2 = Lambda_rp*(1-(rr2*a2))
    return  Q1/Q2

def threshold (X_c1,tau_ph1):
    xc= 5*X_c1
    norma_f1 = (2*gamma*vg*tau_ph1)**(-0.5)
    a = -1*xc/mu
    b = (1 + (r**2))**0.5 - (2*r)
    c = -1 
    delta = (b**2) - (4*a*c)
    ap = ((-1*b) + np.sqrt(delta))/(2*a)
    return (((ap*norma_f1**2)*Aeff))

Delta = np.linspace(0, 10, 1000)
cou = np.linspace(10, 300, 1000)
cou = cou*1e-12
Q = Qfactor(Radius*1e4,cou,1e-3)

y1 =  np.linspace(0, 10, 1000)
y2 =  np.linspace(0, 10, 1000)
y3 =  np.linspace(0, 10, 1000)
y4 =  np.linspace(0, 10, 1000)
y5 =  np.linspace(0, 10, 1000)

for i in range(0,len(cou)):
    y1[i] = threshold (0.01,cou[i])
    y2[i] = threshold (0.05,cou[i])
    y3[i] = threshold (0.1,cou[i])
    y4[i] = threshold (0.2,cou[i])
    y5[i] = threshold (0.4,cou[i])

cou = cou*1e12
plt.plot(cou,y1)
plt.plot(cou,y2,label= r"Radius R = 50 $\mu$ m")
plt.plot(cou,y3,label= r"Radius R = 75 $\mu$ m")
plt.plot(cou,y4,label= r"Radius R = 100 $\mu$ m")
plt.plot(cou,y5,label= r"Radius R = 250 $\mu$ m")
plt.legend()
plt.ylabel(r"Circulating power threshold $P_{th}$ (W)")
plt.xlabel(r"Photon lifetime $\tau_{ph} (ps)$")
plt.xscale("log")






