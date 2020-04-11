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

from time import time

def H(a,x):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM a
    and Gaussian component 2 sigma**2 = one.
    """
    return np.real(wofz(x + 1j*a))


sigma_0 = np.pi * const.e.gauss**2 / const.m_e /const.c


def R(lam1,N,b,line):

    lam_0 = line.lam_0 * unit.Angstrom
    om_0 = 2 * np.pi * const.c / lam_0
    f = line.f
    gamma = line.gamma * unit.s**-1
    
    N = N * unit.cm**-2
    b = b * unit.km/unit.s
    
    lam = lam1 * unit.Angstrom + lam_0
    
    a =  np.float((gamma / ( 4 * om_0) * const.c / b).decompose()) 
    
    om = 2 * np.pi * const.c / lam
    u = np.array(((om - om_0)/ om_0 * const.c / b).decompose())
    
    sigma = 2 * np.sqrt(np.pi) * lam * f * sigma_0 / b * H(a,u)
    
    tau = sigma * N
        
    return 1- np.exp(-tau)

def eq_width(N,b,line):
    return integrate.quad(lambda lam1 : R(lam1,N,b,line), -np.inf,np.inf)


t1 = time()
W = eq_width(12e14,20,lines[0])
t2 = time()

#N = 12e14
#b = 20
#
#line = lines[0]
#
#
#lam_0 = line.lam_0 * unit.Angstrom
#om_0 = 2 * np.pi * const.c / lam_0
#f = line.f
#gamma = line.gamma * unit.s**-1
#
#N = N * unit.cm**-2
#b = b * unit.km/unit.s
#
#a =  np.float((gamma / ( 4 * om_0) * const.c / b).decompose()) 
#
#om = 2 * np.pi * const.c / lam
#u = np.array(((om - om_0)/ om_0 * const.c / b).decompose())
#
#sigma = 2 * np.sqrt(np.pi) * lam * f * sigma_0 / b * H(a,u)
#
#tau = sigma * N


def eq_width_trapz(N,b,line):
    lam = np.linspace(1200.,1230,4000)
    lam1 = lam - line.lam_0
    
    R_list = R(lam1,N,b,line)
    return float(integrate.trapz(R_list, lam))
    

t3 = time()
W2 = eq_width_trapz(12e14,20,lines[0])
t4 = time()


print(t2-t1,t4-t3)






















