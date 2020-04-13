# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 01:04:32 2020

@author: premv
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as unit
import astropy.constants as const

from voigt_and_W import eq_width, W_by_lam_linear, W_by_lam_saturated, W_by_lam_damped, plot_CoG

c = float(const.c/(unit.m/unit.s))

from analyse_spectrum import lines 



line = lines['HI'][0]




#sigma_0_dim = np.pi * const.e.gauss**2 / const.m_e /const.c
#sigma_0 = np.float(sigma_0_dim.decompose()/unit.m**2*unit.s)

    


#def plot_CoG(N_list,b_list,line,gen_plot=True,save_file=None,check_approx=True):
#    plt.figure(dpi=200,figsize=(9,7))
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
#        if check_approx:
#            N_lam_f = N_list * line.lam_0 * line.f
#            plt.plot(N_lam_f[5:16], W_by_lam_saturated(N_lam_f[5:16],b),'--',label="saturated b = {} km/s".format(b))
#            
#    if check_approx:
#        N_lam_f = N_list * line.lam_0 * line.f
#        plt.plot(N_lam_f[:10], W_by_lam_linear(N_lam_f[:10]),'--',label="linear approx")
##        plt.plot(N_lam_f, W_by_lam_saturated(N_lam_f,b),'--',label="saturated")
#        plt.plot(N_lam_f[6:], W_by_lam_damped(N_lam_f[6:],line),'--',label="damped approx")
#    plt.ylim(2e-7,1e-1)
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
#    
#    
#    plt.close()


line_list = lines['HI']

N_list = np.logspace(11,22,23) #* unit.cm**-2
#NH = NH_list[9]
#NH = 1 * 10**12.5 * unit.cm**-2
#b_list = np.logspace(0,2,10) #* unit.km/unit.s
b_list = 2**np.linspace(0,6,7) #* unit.km/unit.s
#b = b_list[1]

#plot_CoG(N_list,b_list,line,save_file='CoG_Ly-alpha')
plot_CoG(N_list,b_list,line_list,gen_plot=True,check_approx=True,save_file='CoG_Ly-alpha-approx')


#W_by_lam_appr = N_list[:5] * line.lam_0 * line.f / 1.132e20
#W_appr = N_list[:5] * line.lam_0**2 * line.f / 1.132e20


#plt.plot(N_list[:5], W_appr)

#factor = NH_list * line.lam_0 * line.f / (W_by_lam)



#N_lam_f = 10**np.linspace(14,25,120)
#
#plt.plot(N_lam_f, W_by_lam_linear(N_lam_f),'--',label="linear")
#plt.plot(N_lam_f, W_by_lam_saturated(N_lam_f,16),label="saturated")
#plt.plot(N_lam_f, W_by_lam_damped(N_lam_f,line),label="damped")
#
#plt.legend()
#plt.xscale('log')
#plt.yscale('log')










