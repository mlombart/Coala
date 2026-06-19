import sympy as sym
import numpy as np
from matplotlib import pyplot as plt

import os
import re
import sys

sys.path.append("../..")

from functions_plot import *



#options for plots
plt.rcParams["font.size"]= 14
plt.rcParams['lines.linewidth'] = 3
plt.rcParams["legend.columnspacing"] = 0.5

marker_style = dict( marker='o',markersize=8, markerfacecolor='white', linestyle='',markeredgewidth=2)
savefig_options=dict(bbox_inches='tight')

########################################################
#/!\ need to be the same value used in coala_hydro.f90 #
########################################################
nbins = 50
dtg = 1e-2
rho_gas = 1e-15
coeff_pl = -3.5
##################################################

source_dv = 'dv_ormel+dv_brownian'
path_data = '../data/' + source_dv


#load data

list_kpol = get_order_kpol(path_data+'/nbins='+str(nbins))

dict_data = {}
print("Loading data for order k =",(list_kpol))
for k in list_kpol:
   dict_data['order k=%d'%(k)] = {}

   dict_data['order k=%d'%(k)]['kpol']          = k
   dict_data['order k=%d'%(k)]['sizegrid']      = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/sizegrid.txt')
   dict_data['order k=%d'%(k)]['sizemeanlog']   = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/sizemeanlog.txt')
   dict_data['order k=%d'%(k)]['time']          = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/time.txt')
   dict_data['order k=%d'%(k)]['rhodust_t0']    = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/rhodust_t0.txt')
   dict_data['order k=%d'%(k)]['rhodust_tend']  = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/rhodust_tend.txt')
   

   #compute mass fraction

   rhodust_t0   = dict_data['order k=%d'%(k)]['rhodust_t0'] 
   rhodust_tend = dict_data['order k=%d'%(k)]['rhodust_tend'] 

   mass_frac_t0 = rhodust_t0/np.sum(rhodust_t0)
   mass_frac_tend = rhodust_tend/np.sum(rhodust_tend)

   dict_data['order k=%d'%(k)]['mass_fraction_t0']   = mass_frac_t0
   dict_data['order k=%d'%(k)]['mass_fraction_tend'] = mass_frac_tend




#define limits for plot
first_key = list(dict_data.keys())[0]
smin = dict_data[first_key]['sizegrid'][0]
smax = dict_data[first_key]['sizegrid'][-1]
tend = dict_data[first_key]['time'][-1]

ymin = 10**(-16)
ymax = 10

yr = 31556926. #cgs


fig,axs = plt.subplots(1,2,figsize=(12,6))

for i in dict_data:
   
   kpol              = dict_data[i]['kpol']
   sizemeanlog       = dict_data[i]['sizemeanlog']
   mass_frac_t0      = dict_data[i]['mass_fraction_t0']
   mass_frac_tend    = dict_data[i]['mass_fraction_tend']
   rhodust_t0        = dict_data[i]['rhodust_t0']
   rhodust_tend      = dict_data[i]['rhodust_tend']
   tff           = dict_data[i]['time'][0]


   if (kpol==0):
      color = 'black'
      
   else:
      color = "C"+str(kpol)

   cm_to_mum = 1e4

   total_dust_density_t0   = np.sum(rhodust_t0)
   total_dust_density_tend = np.sum(rhodust_tend)
   abs_err = np.abs(total_dust_density_tend-total_dust_density_t0)/total_dust_density_t0

   #plot rhodust
   axs[0].loglog(sizemeanlog*cm_to_mum,rhodust_t0,markeredgecolor=color,**marker_style,alpha=0.5,zorder=4)
   axs[0].loglog(sizemeanlog*cm_to_mum,rhodust_tend,markeredgecolor=color,label='order k=%d, abs error mass = %.2e'%(kpol,abs_err),**marker_style,zorder=4)

   #plot mass fraction
   axs[1].loglog(sizemeanlog*cm_to_mum,mass_frac_t0,markeredgecolor=color,**marker_style,alpha=0.5,zorder=4)
   axs[1].loglog(sizemeanlog*cm_to_mum,mass_frac_tend,markeredgecolor=color,**marker_style,zorder=4)


axs[0].set_xlim(smin*cm_to_mum,smax*cm_to_mum)
# plt.ylim(1e-30,1e-10)
axs[0].set_xlabel(r'size $\mu m$')
axs[0].set_ylabel(r'dust density')

axs[0].legend(loc='lower left',ncol=1)


axs[1].set_xlim(smin*cm_to_mum,smax*cm_to_mum)
axs[1].set_ylim(1e-7,1.)
axs[1].set_xlabel(r'size $\mu m$')
axs[1].set_ylabel(r'mass fraction')
axs[1].yaxis.tick_right()
axs[1].yaxis.set_label_position("right")



fig.suptitle('source ' +source_dv + r', $\tau$=%d tff'%(tend/tff))
plt.tight_layout()


# to save plot
# plt.savefig('./source_' + source_dv + '.pdf',**savefig_options)

plt.show()

