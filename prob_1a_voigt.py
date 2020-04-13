# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 23:32:06 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
#import read_atoms
from voigt_and_W import get_voigt
from analyse_spectrum import lines 



#lines = read_atoms.read_lines('file/atom_identified.dat')
line = lines['HI'][0]



#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s

#N = N_list[9]
N = 1 * 10**13.5 #* unit.cm**-2

plt.figure(dpi=150,figsize=(7,5))

for b in b_list:
    voigt = get_voigt(N,b,line,1215,1216.5)
    plt.plot(voigt[0],1-voigt[1],label="{} km/s".format(b))
    
plt.title(r"Voigt Profile for Ly-$\alpha$ with N=10$^{%.1f}$cm$^{-2}$"%np.log10(N))
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"Normalised flux")
plt.legend()
plt.savefig('Voigt-Ly-a_vary_b.pdf',bbox_inches='tight')
  
plt.close()
    
N_list = np.logspace(13,20,8) #* unit.cm**-2
b = 12
    
plt.figure(dpi=150,figsize=(7,5))    
for N in N_list:
    voigt = get_voigt(N,b,line,1215,1216.5)
    plt.plot(voigt[0],1-voigt[1],label="N=10$^{%.0f}$cm$^{-2}$"%np.log10(N))
    
plt.title(r"Voigt Profile for Ly-$\alpha$ with b={} km/s".format(b))
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"Normalised flux")
plt.legend()

plt.savefig('Voigt-Ly-a_vary_N.pdf',bbox_inches='tight')


