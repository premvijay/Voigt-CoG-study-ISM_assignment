# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 01:05:15 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt
#import read_atoms
from voigt_and_W import eq_width, plot_CoG
from analyse_spectrum import lines



#def plot_CoG(N_list,b_list,line_list,gen_plot=True,save_file=None,check_approx=False):
#    plt.figure(dpi=200,figsize=(9,7))
#    for line in line_list:
#        for b in b_list:
#            W = np.zeros(len(N_list))
#            for i in range(len(N_list)):
#                W[i] = (eq_width(N_list[i],b,line))
#        #    W = np.array(W) 
#            if gen_plot:
#                plt.plot(N_list*line.f*line.lam_0, W/line.lam_0, label="{}-{:.0f}-line with b={} km/s".format(line.ID,line.lam_0,b))
#            else:
#                plt.plot(N_list, W, label="{}-{:.0f}-line with b={} km/s".format(line.ID,line.lam_0,b))
#    #            pass
#            if check_approx:
#                N_lam_f = N_list * line.lam_0 * line.f
#                plt.plot(N_lam_f[5:16], W_by_lam_saturated(N_lam_f[5:16],b),'--',label="saturated b = {} km/s".format(b))
#        if check_approx:
#            plt.plot(N_lam_f[6:], W_by_lam_damped(N_lam_f[6:],line),'--',label="damped approx")
#
#    if check_approx:
##        N_lam_f = N_list * line.lam_0 * line.f
#        plt.plot(N_lam_f[:10], W_by_lam_linear(N_lam_f[:10]),'--',label="linear approx")  
#
#    plt.legend(framealpha=0.2)    
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
##    plt.close()
#    return plt.gcf(), plt.gca()
    


line_list = lines['HI'] + lines['SiII']

N_list = np.logspace(11,22,23) #* unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 10**np.linspace(0,1,2) #* unit.km/unit.s
#b = b_list[1]

plot_CoG(N_list,b_list,line_list,save_file='CoG_multi',legend_id=True)
plot_CoG(N_list,b_list,line_list,gen_plot=False,save_file='CoG_multi',legend_id=True)