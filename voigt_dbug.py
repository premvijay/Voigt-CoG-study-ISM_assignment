# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:33:18 2020

@author: premv
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import astropy.units as unit
from astropy.constants import c as c1#, e, m_e #const
#from astropy.modeling.models import Voigt1D

from scipy.special import wofz

def Voigt2(a,x):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM a
    and Gaussian component 2 sigma**2 = one.
    """
    return np.real(wofz(x + 1j*a))



c = 3e8
pi = np.pi

@np.vectorize
def Voigt(a, u):
    I = integrate.quad(lambda y: np.exp(-y**2)/(a**2 + (u - y)**2),-np.inf, np.inf)[0]

    return (a/pi)*I

from time import time

#sigma_0 = pi * e.gauss**2 / m_e /c

#NH_list = np.logspace(11,22,12) * unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
b_list = np.linspace(10,50,5) * 1000#unit.km/unit.s
b_list1 = np.linspace(10,50,5) * unit.km/unit.s
b = b_list[1]
b1 = b_list1[1]

lam_0 = 1215.6701e-10# * unit.Angstrom
om_0 = 2 * pi * c / lam_0
f = 0.4164
gamma = 6.265 * 10**8 #* unit.s**-1


lam_01 = 1215.6701 * unit.Angstrom
om_01 = 2 * pi * c1 / lam_01
f = 0.4164
gamma1 = 6.265 * 10**8 * unit.s**-1



#Del_om_D = 2 * pi * b / lam_0
a =  np.float((gamma / ( 4 * om_0) * c / b))#.decompose()) 
a1 =  np.float((gamma1 / ( 4 * om_01) * c1 / b1).decompose())


def H1(a,u):
    return Voigt2(a,u)

#H2 = Voigt1D(x_0=0, amplitude_L=.57/a, fwhm_L=2*a, fwhm_G=2*np.sqrt(np.log(2)))


lam = np.linspace(1215,1216,400) * 1e-10# unit.Angstrom
om = 2 * pi * c / lam

lam1 = np.linspace(1215.665,1215.675,400) * unit.Angstrom
om1 = 2 * pi * c1 / lam1


u = np.array((om - om_0)/ om_0 * c / b)
u1 = np.array((om1 - om_01)/ om_01 * c1 / b1)

#sigma = 2 * np.sqrt(pi) * lam * f * sigma_0 / b * H1(a,u)

phi = 2 * np.sqrt(pi) * lam / b * H1(a,u)
phi1 = 2 * np.sqrt(pi) * lam1 / b1 * H1(a,u)

normalize = integrate.trapz(phi, om)/ (2*pi)
print(normalize)

#plt.plot(om,phi)



#Gamma = 2 * e.gauss**2 / (3 * m_e * c)  * 4 * pi**2 * lines[3].f / (lines[3].lam_0 * unit.Angstrom)**2 


#a = 1
#u = np.linspace(-30,30,20000)
#
#t1 = time()
#V = Voigt(a,u)
#t2 = time()
#
#V2 = Voigt2(a,u)
#t3 = time()
#
#plt.plot(u,V)
#plt.plot(u,V2)
#
#I = integrate.trapz(V2,u)
#
#print(t2-t1,t3-t2)










