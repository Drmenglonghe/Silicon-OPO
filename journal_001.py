import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from scipy import integrate
import matplotlib.pyplot as plt

##########################################################################################################################
def Linear_response(r1,alpha_dB1):
    alpha1 = np.log(10)*alpha_dB1/10
    return alpha1 

def photonlifetime(r2,theta2,alpha_dB1):
    r2 = r2
    tR   = (2 * np.pi * r2) / vg
    L = 2*np.pi*r2
    Loss = Linear_response(r2,alpha_dB1)*L + theta2
    return  tR/Loss

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

##########################################################################################################################
X_c1 = 4
r = 0.2
mu = 25
Aeff = 1e-9
#delta_T 
sigmma = 1.45e-17 #FCA corss section 
Beta_GVD = -0.62*(1e-24/1e2)#, 0.1*(1e-24/1e2)]

##########################################################################################################################
c = 3e10 # cm/s
ng = 3.6 
vg = c/ng 
X_th = 0.34
Tau_th = 0.04 * 1e-6  
Aeff = 0.1e-8 
neff = Aeff
theta = 1e-1# critical coupling 
radius = 50e-4 # um in cm
alpha_dB = 1 # dB/cm
gamma = 3.1e-9

tR = (2 * np.pi * radius) / vg # s 
tau_ph = photonlifetime(radius,theta,alpha_dB) # s 

normal_f = (2*gamma*vg*tau_ph)**(-0.5)
normal_in = ((tR**2)/(8*theta*gamma*vg*(tau_ph**3)))**(0.5)
normal_t = (abs(Beta_GVD)*vg*tau_ph)**(0.5)
normal_tau = 2*tau_ph  
normal_n = 1/(mu*sigmma*vg*tau_ph)
#normal_T = 1/(2*delta_T*tau_ph*vg)

##########################################################################################################################
f0 = (1400e-3)/(normal_f**2*neff)
Field1 = np.linspace(0,f0,1000)
d0 = 1000e9*(2*np.pi)*(2*tau_ph)
delta1 = np.linspace(-1*d0,1*d0,1000)

delta,Field=np.meshgrid(delta1,Field1)
z1 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (10e-3)/(normal_in**2*neff)
z2 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (20e-3)/(normal_in**2*neff)
z3 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
 )**2))  - (30e-3)/(normal_in**2*neff)
z4 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (40e-3)/(normal_in**2*neff)
z5 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (50e-3)/(normal_in**2*neff)
z6 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (60e-3)/(normal_in**2*neff)
z7 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (70e-3)/(normal_in**2*neff)
z8 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (80e-3)/(normal_in**2*neff)
z9 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (90e-3)/(normal_in**2*neff)
z10 = Field * (((1+ (r*Field) + (X_c1*(Field**2)/mu))**2) +((delta-  (Field) +  (X_c1*(Field**2)) - ((X_th*Tau_th)/(2*tau_ph))*((r*(Field**2)) + ((2*X_c1*(Field**3)/mu)))
)**2))  - (100e-3)/(normal_in**2*neff)

for i in range(0,len(Field1)):
    for j in range(0, len(delta1)):
        delta[i][j] = delta[i][j]/(1e9*(2*np.pi)*(2*tau_ph))
                
for i in range(0,len(Field1)):
    for j in range(0, len(delta1)):
        Field[i][j] = ((Field[i][j]*(normal_f**2)*Aeff))

plt.contour(delta,Field,z1,0,colors=lighten_color('r', 0.2))
#plt.contour(delta,Field,z2,0,colors=lighten_color('r', 0.2))
plt.contour(delta,Field,z3,0,colors=lighten_color('r', 0.4))
#plt.contour(delta,Field,z4,0,colors=lighten_color('r', 0.4))
plt.contour(delta,Field,z5,0,colors=lighten_color('r', 0.6))
#plt.contour(delta,Field,z6,0,colors=lighten_color('r', 0.6))
plt.contour(delta,Field,z7,0,colors=lighten_color('r', 0.8))
#plt.contour(delta,Field,z8,0,colors=lighten_color('r', 0.8))
plt.contour(delta,Field,z9,0,colors=lighten_color('r', 1))
#plt.contour(delta,Field,z10,0,colors=lighten_color('r', 1))

plt.xlabel(r"Cold detuning ${\Delta}_{0}$ (GHz)")
plt.ylabel(r"Circulating field power $P_{c}$ (W)")
plt.savefig("journal14.png",dpi=1000,bbox_inches='tight')
plt.show
