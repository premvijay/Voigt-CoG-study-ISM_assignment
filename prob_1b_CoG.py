# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 01:04:32 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
#import read_atoms
from voigt_and_W import eq_width, plot_CoG
from analyse_spectrum import lines 



#plt.close()


#def plot_CoG(N_list,b_list,line,gen_plot=True,save_file=None):
#    plt.figure(dpi=150)
#    for b in b_list:
#        W = np.zeros(len(N_list))
#        for i in range(len(N_list)):
#            W[i] = (eq_width(N_list[i],b,line))
#    #    W = np.array(W) 
#        if gen_plot:
#            plt.plot(N_list*line.f*line.lam_0, W/line.lam_0, label="b={} km/s".format(b))
#        else:
#            plt.plot(N_list, W, label="b={} km/s".format(b))
##            pass
#
#    
#    plt.legend()    
#    plt.xscale('log')
#    plt.yscale('log')
#    if gen_plot:
#        plt.xlabel(r'$ N_{H} \lambda f (cm^{-2} \AA)$')
#        plt.ylabel(r'$ W_{\lambda} / \lambda $')
#        plt.title("Curve of Growth")
#        if save_file is not None:
#            plt.savefig(save_file+"_gen.pdf",bbox_inches='tight')
#    else:
#        plt.xlabel(r'$ N_{H} (cm^{-2})$')
#        plt.ylabel(r'$ W_{\lambda} (\AA)$')
#        plt.title("Curve of Growth")
#        if save_file is not None:
#            plt.savefig(save_file+".pdf",bbox_inches='tight')
#    #plt.figure(figsize=(13,8))
#    plt.close()


line_list = lines['HI']

N_list = np.logspace(11,22,23) #* unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s
#b = b_list[1]

plot_CoG(N_list,b_list,line_list,save_file='CoG_Ly-alpha')
plot_CoG(N_list,b_list,line_list,gen_plot=False,save_file='CoG_Ly-alpha')

 
#    if gen_plot:
#W_by_lam_appr = N_list[:5] * line.lam_0 * line.f / 1.132e20
#W_appr = N_list[:5] * line.lam_0**2 * line.f / 1.132e20

#plt.plot(N_list[:5] * line.lam_0 * line.f, W_by_lam_appr)
#plt.plot(N_list[:5], W_appr)

#factor = NH_list * line.lam_0 * line.f / (W_by_lam)
    
    
