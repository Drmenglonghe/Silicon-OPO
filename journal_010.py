# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:21:56 2022

@author: mengl
"""

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
from mpl_toolkits.axes_grid1 import make_axes_locatable

X_c = 2
r = 0.189 
mu = 25 

c = 3e10                             # Speed of light in cm / s
gamma = 3.1e-9                       # cm/W
ng = 3.6                             # Refractive index of Silicon at 1.55um. 
vg = c/ng                            # Group velocity
theta = 0.1                         # Assuming a set % of power is coupled to ring from bus 
Lambda_rp = 1550e-7                  # cm 1550 nm 
wr = 2 * np.pi * c / Lambda_rp       # rad/s 1550 nm
alpha_dB = 1 
Radius = 100e-4                       # good choose for critical coupling 
tR = (2 * np.pi * Radius) / vg
Aeff = 1e-9
wr = 2 * np.pi * c / Lambda_rp  
fsr = Lambda_rp**2/(ng*2*np.pi*Radius)*1e7

def Linear_response(r1,alpha_dB1):
    alpha1 = np.log(10)/10*alpha_dB1
    return alpha1

def photonlifetime(r2,theta2,alpha_dB1):
    r2 = r2
    tR   = (2 * np.pi * r2) / vg
    L = 2*np.pi*r2
    Loss = Linear_response(r2,alpha_dB1)*L + theta2
    return  tR/Loss

def Qfactor(r2,alpha,theta2):
    r2 = r2
    rr2 = (1-theta2)**0.5
    L = 2*np.pi*r2
    a2 = np.exp(-1*Linear_response(r2,alpha)*L)**0.5
    L2 = 2*np.pi*r2
    Q1 = np.pi * ng * L2*np.sqrt(a2*rr2)
    Q2 = Lambda_rp*(1-(rr2*a2))
    return  Q1/Q2

def threshold (tau_ph1,X_c1):
    tau_ph1 = tau_ph1*1e-12
    X_c1 = 5* X_c1
    xc = X_c1
    norma_f1 = (2*gamma*vg*tau_ph1)**(-0.5)
    a = -1*xc/mu
    b = (1 + (r**2))**0.5 - (2*r)
    c = -1 
    delta = (b**2) - (4*a*c)
    ap = ((-1*b) + np.sqrt(delta))/(2*a)
    return (((ap*norma_f1**2)*Aeff))

x = np.linspace(50, 500, 100)
y = np.linspace(0.05, 0.5, 100)
Z0 = np.zeros((len(y),(len(x))))

for i in range(0,len(y)):
    for j in range(0,len(x)):
        Z0[i][j] = threshold(x[j],y[i])

fig, ax = plt.subplots()
im = ax.contourf(x, y, Z0, 24,cmap = "hot_r")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Photon lifetime $\tau_{ph}$ (ps)")
plt.ylabel(r"Free carrier/photon lifetime ($\tau{c}$/$\tau{ph}$)")
plt.title(r" Circulating power threshold $P_{th}$ (W)")
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size = '4%', pad = 0.5)
fig.add_axes(cax)
cbar = fig.colorbar(im, cax = cax, orientation = 'horizontal')
plt.savefig("dopo001.png",dpi=1000,bbox_inches='tight')
plt.show()






