# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 17:23:15 2020

@author: premv
"""

import numpy as np


from analyse_spectrum import lines
#from voigt_and_W import eq_width
from voigt_and_W import plot_CoG

lines_Fe = lines['FeII']

N_list = np.logspace(13,19,13) #* unit.cm**-2
#NH = NH_list[9]
#b_list = 2**np.linspace(0,5,6) #* unit.km/unit.s
b_list = np.linspace(18,24,1)
#b = b_list[1]

fig,ax = plot_CoG(N_list,b_list,lines_Fe[:1],legend_id=True)
#plot_CoG(N_list,b_list,line_list,gen_plot=False,save_file='CoG_multi',legend_id=True)

N_list = np.logspace(15.25,15.3,1)

for N in N_list:
    for line in lines_Fe:
        ax.scatter(N*line.f*line.lam_0,line.W_by_lam, label="{}".format(N))



lines_Ni = lines['NiII']

N_list = np.linspace(1.4,1.5,1) * 1e14

for N in N_list:
    for line in lines_Ni:
        ax.scatter(N*line.f*line.lam_0,line.W_by_lam, label="{}".format(N))

fig.savefig("FeII-NiII-consist-fiting.pdf",bbox_inches='tight')




