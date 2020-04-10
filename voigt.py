# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 19:16:33 2020

@author: premv
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import astropy.units as unit
import astropy.constants as const
from astropy.modeling.models import Voigt1D

from time import time

@np.vectorize
def Voigt(a, u):
    I = integrate.quad(lambda y: np.exp(-y**2)/(a**2 + (u - y)**2),-np.inf, np.inf)[0]

    return (a/np.pi)*I

#def VoigtRange(a, u):
#    Integ = 
#    I = integrate.quad(lambda y: np.exp(-y**2)/(a**2 + (u - y)**2),-np.inf, np.inf)[0]
#
#    return (a/np.pi)*I

a = 1
u = np.linspace(-10,10,20)

plt.figure()

t1 = time()
V2 = Voigt(a,u)


t2 = time()
plt.plot(u,V2)

#x = u#np.arange(0, 10, 0.01)
v1 = Voigt1D(x_0=0, amplitude_L=.57/a, fwhm_L=2*a, fwhm_G=2*np.sqrt(np.log(2)))

t3 = time()
V1 = v1(u)

t4 = time()

plt.plot(u, V1,label="astropy")

plt.legend()
plt.show()

dV = V2-V1
fV = V2/V1

print(t2-t1,t4-t3)

plt.close()

#class Voigt:
#    def __init__(self):
#        self.phi = 0
        

sigma_0 = np.pi * const.e.gauss**2 / const.m_e /const.c

NH_list = np.logspace(11,22,12) * unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
b_list = np.linspace(10,50,5) * unit.km/unit.s
#b = b_list[1]

lam_0 = 1215.6701 * unit.Angstrom
om_0 = 2 * np.pi * const.c / lam_0
f = 0.4164
gamma = 6.265 * 10**8 * unit.s**-1 #* np.pi


for b in b_list:
    W = []
    Del_om_D = 2 * np.pi * b / lam_0
    a =  np.float((gamma / ( 4 * om_0) * const.c / b).decompose()) 
    
    
    
    def H1(a,u):
        return Voigt(a,u)
    
    H2 = Voigt1D(x_0=0, amplitude_L=.57/a, fwhm_L=2*a, fwhm_G=2*np.sqrt(np.log(2)))
        
    
    lam = np.linspace(1210.,1220,400) * unit.Angstrom
    om = 2 * np.pi * const.c / lam
    
    u = np.array(((om - om_0)/ om_0 * const.c / b).decompose())
    
    sigma = 2 * np.sqrt(np.pi) * lam * f * sigma_0 / b * H1(a,u)
    
    for NH in NH_list:
        tau = sigma * NH
        
        flux = np.exp(-tau)
        
#        plt.plot(lam,flux,label="b={},N={:.3e}".format(b,NH))
        #plt.plot(u,H2(u))
        
        W.append(float(integrate.trapz(1 - flux, lam)/unit.Angstrom))
    
    
#    plt.legend()
    
    W = np.array(W)
#    plt.savefig("multi_voigt.pdf")
#    plt.close()
    
    #plt.plot(NH_list,W)
    #
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel(r'$\log [ N_{H} (cm^{-2})]$')
    #plt.ylabel(r'$\log [ W_{\lambda} (\AA)]$')
    #
    plt.plot(NH_list*f*lam_0, W/lam_0*unit.Angstrom, label="b={}".format(b))
    

plt.legend()    
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\log [ N_{H} \lambda f (cm^{-2} \AA)]$')
plt.ylabel(r'$\log [ W_{\lambda} / \lambda]$')











