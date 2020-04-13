# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:13:28 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as unit
import astropy.constants as const


from analyse_spectrum import lines, spectrum
from voigt_and_W import get_voigt

line = lines['HI'][0]

c = float(const.c/(unit.m/unit.s))


sigma_0_dim = np.pi * const.e.gauss**2 / const.m_e /const.c
sigma_0 = np.float(sigma_0_dim.decompose()/unit.m**2*unit.s)

plt.figure(figsize=(9,7))

plt.plot(spectrum[:,0],spectrum[:,3]/1.6e-14,label=r"Smoothened spectrum Ly-$\alpha$")
plt.xlim(1150,1270)
#ax.set_ylim(0,1.1)

linear_factor = (np.pi* const.e.gauss**2/const.m_e / const.c**2).decompose()

#def W_by_lam_damped(N_lam_f,line):
#    return 2*np.sqrt(N_lam_f * linear_factor * line.gamma * line.lam_0 * 1e-10 /c/np.pi**3)

#def N_lam_f_damped(W_by_lam,line):

W = line.W * unit.Angstrom
lam = line.lam_0 * unit.Angstrom
gamma = line.gamma * unit.s**-1
f = line.f


line.set_N(float(W**2 * np.pi**3 * const.c / (4 * lam**4 * f * linear_factor * gamma) /unit.cm**-2))


print(linear_factor,line.N)

b_list = 10**np.linspace(0,2,3)

#plt.figure(dpi=150,figsize=(7,5))

for b in b_list:
    voigt = get_voigt(line.N,b,line,1150,1270)
    plt.plot(voigt[0],1-voigt[1],label="Voigt profile with b = {} km/s".format(b),linestyle='--')
    
plt.title(r"Voigt Profile overlay on spectrum for Ly-$\alpha$ with N=%.3e cm$^{-2}$"%(line.N))
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"Normalised flux")
plt.legend()



#plt.scatter((line.W),(line.N))

plt.savefig('Voigt-Ly-a_overlay.pdf',bbox_inches='tight')







