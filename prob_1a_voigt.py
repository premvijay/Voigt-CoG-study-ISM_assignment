# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 23:32:06 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt



N_list = np.logspace(13,20,8) #* unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s

line = lines[0]

#N = N_list[9]
N = 1 * 10**13.5 #* unit.cm**-2

plt.figure(dpi=150)

for b in b_list:
    voigt = get_voigt(N,b,line,1215,1217)
    plt.plot(voigt[0],1-voigt[1])
  
plt.close()
    
b = b_list[3]
    
plt.figure(dpi=150)    
for N in N_list:
    voigt = get_voigt(N,b,line,1215,1217)
    plt.plot(voigt[0],1-voigt[1])



