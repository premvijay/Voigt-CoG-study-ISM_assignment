# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 01:04:32 2020

@author: premv
"""

import numpy as np
import read_atoms

from voigt_and_W import plot_CoG

lines = read_atoms.read_lines('file/atom_identified.dat')
line_list = lines['HI']

N_list = np.logspace(11,22,23) #* unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s
#b = b_list[1]

#plot_CoG(N_list,b_list,line,save_file='CoG_Ly-alpha')
plot_CoG(N_list,b_list,line_list,gen_plot=True,check_approx=True,save_file='CoG_Ly-alpha-approx',fig_size=(9,7))









