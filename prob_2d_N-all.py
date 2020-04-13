# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:28:35 2020

@author: premv
"""

#import numpy as np



from analyse_spectrum import lines
#from voigt_and_W import eq_width
from voigt_and_W import linear_factor


for ID in lines:
    if ID != 'HI' and ID != 'FeII' and ID != 'NiII':
        N_mean = 0
        for line in lines[ID]:
            line.set_N(line.W_by_lam / line.lam_0 / line.f / linear_factor)
    #        print("ID {} W_by_lam = {:.3e} and N = {:.3e} cm**-2".format(line.ID,line.W_by_lam,line.N))
            N_mean += line.N
        N_mean /= len(lines[ID])
        print("ID {} W_by_lam = {:.3e} and N = {:.3e} cm**-2".format(line.ID,line.W_by_lam,line.N))
        
        





