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

#from time import time

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
    return integrate.quad(R, -np.inf,np.inf,args=(N,b,line))[0]


def get_voigt(N,b,line,start,stop):
    lam = np.linspace(start,stop,4000)
    lam1 = lam - line.lam_0
    
    return lam,R(lam1,N,b,line)

def eq_width_trapz(N,b,line,start,stop):    
    voigt = get_voigt(N,b,line,start,stop)
#    plt.plot(lam,1-R_list)
    return float(integrate.trapz(1-voigt[1], voigt[0]))
    

#N_list = np.logspace(11,22,23) #* unit.cm**-2
##NH = NH_list[9]
##NH = 1 * 10**12.5 * unit.cm**-2
##b_list = np.logspace(0,2,10) #* unit.km/unit.s
#b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s
##b = b_list[1]
#
#from analyse_spectrum import lines
#line = lines['HI'][0]

#for b in b_list:
#    W = []
#    for N in N_list:
#        W.append(eq_width(N,b,line))
#    W = np.array(W) 
##    plt.plot(N_list*line.f*line.lam_0, W/line.lam_0, label="b={} km/s".format(b))
#    plt.plot(N_list, W, label="b={} km/s".format(b))









