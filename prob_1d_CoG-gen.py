# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 01:05:15 2020

@author: premv
"""

import numpy as np
import read_atoms
from voigt_and_W import plot_CoG

    

lines = read_atoms.read_lines('file/atom_identified.dat')
line_list = lines['HI'] + lines['SiII']

N_list = np.logspace(11,22,23) #* unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 10**np.linspace(0,1,2) #* unit.km/unit.s
#b = b_list[1]

plot_CoG(N_list,b_list,line_list,save_file='CoG_multi',legend_id=True)
plot_CoG(N_list,b_list,line_list,gen_plot=False,save_file='CoG_multi',legend_id=True)