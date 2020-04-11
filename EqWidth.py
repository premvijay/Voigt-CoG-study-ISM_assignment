# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 13:04:16 2020

@author: premv
"""

import numpy as np
from scipy import integrate
from scipy.special import wofz
import matplotlib.pyplot as plt
import astropy.units as unit
import astropy.constants as const

c = float(const.c/(unit.m/unit.s))

from time import time

def H(a,x):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM a
    and Gaussian component 2 * sigma**2 = one.
    """
    return np.real(wofz(x + 1j*a))


sigma_0_dim = np.pi * const.e.gauss**2 / const.m_e /const.c
sigma_0 = np.float(sigma_0_dim.decompose()/unit.m**2*unit.s)


def R(lam1,N,b,line):

    lam_0 = line.lam_0 *1e-10#* unit.Angstrom
    om_0 = 2 * np.pi * c / lam_0
    f = line.f
    gamma = line.gamma #* unit.s**-1
    
    N = N * 1e4 #unit.cm**-2
    b = b * 1e3 #unit.km/unit.s
    
    lam = lam1 * 1e-10 + lam_0
    
    a =  np.float((gamma / ( 4 * om_0) * c / b)) 
    
    om = 2 * np.pi * c / lam
    u = np.array(((om - om_0)/ om_0 * c / b))
    
    sigma = 2 * np.sqrt(np.pi) * lam * f * sigma_0 / b / (2*np.pi) * H(a,u)
    
    tau = sigma * N 
        
    return 1-np.exp(-tau)

def eq_width(N,b,line):
    return integrate.quad(1-flux, -np.inf,np.inf,args=(N,b,line))[0]

#
#t1 = time()
#W = eq_width(12e19,20,lines[7])
#t2 = time()

def get_voigt(N,b,line,start=line.start,stop=line.stop):
    lam = np.linspace(start,stop,4000)
    lam1 = lam - line.lam_0
    
    return lam,R(lam1,N,b,line)

def plot_voigt(N,b,line):
    voigt = get_voigt(N,b,line)
    plt.figure(dpi=300)
    plt.plot(voigt[0],1-voigt[1])

def eq_width_trapz(N,b,line):    
    voigt = get_voigt(N,b,line)
#    plt.plot(lam,1-R_list)
    return float(integrate.trapz(1-voigt[1], voigt[0]))
    

t3 = time()
W2 = eq_width_trapz(12e11,200,lines[0])
t4 = time()
plot_voigt(12e12,10,lines[0])
plot_voigt(12e12,20,lines[0])

#print(t2-t1,t4-t3)


N_list = np.logspace(11,22,23) #* unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s
#b = b_list[1]

line = lines[0]

for b in b_list:
    W = []
    for N in N_list:
        W.append(eq_width(N,b,line))
    W = np.array(W) 
#    plt.plot(N_list*line.f*line.lam_0, W/line.lam_0, label="b={} km/s".format(b))
    plt.plot(N_list, W, label="b={} km/s".format(b))



#W_by_lam = W/line.lam_0
#W_by_lam_appr = N_list[:5] * line.lam_0 * line.f / 1.132e20
W_appr = N_list[:5] * line.lam_0**2 * line.f / 1.132e20

#plt.plot(N_list[:5] * line.lam_0 * line.f, W_by_lam_appr)
plt.plot(N_list[:5], W_appr)

#factor = NH_list * line.lam_0 * line.f / (W_by_lam)



plt.legend()    
plt.xscale('log')
plt.yscale('log')
#plt.xlabel(r'$\log [ N_{H} \lambda f (cm^{-2} \AA)]$')
#plt.ylabel(r'$\log [ W_{\lambda} / \lambda]$')
plt.title("Curve of Growth")

plt.xlabel(r'$\log [ N_{H} (cm^{-2})]$')
plt.ylabel(r'$\log [ W_{\lambda} (\AA)]$')

plt.savefig('CoG_Ly-alpha.pdf')
#plt.figure(figsize=(13,8))





















