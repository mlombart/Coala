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



nbins=20
kernel = 'kadd'
path_data = '../data/' + kernel
display_polynomials = True


#load data
# print(os.listdir(path=path_data))
# sys.exit()

list_kpol = get_order_kpol(path_data+'/nbins='+str(nbins))


dict_data = {}
print("Loading data for order k =",(list_kpol))
for k in list_kpol:
   dict_data['order k=%d'%(k)] = {}

   dict_data['order k=%d'%(k)]['kpol']               = k
   dict_data['order k=%d'%(k)]['massgrid']           = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/massgrid.txt')
   dict_data['order k=%d'%(k)]['massbins']           = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/massbins.txt')
   dict_data['order k=%d'%(k)]['massmeanlog']        = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/massmeanlog.txt')
   dict_data['order k=%d'%(k)]['time']               = np.loadtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/time.txt')
   dict_data['order k=%d'%(k)]['gt0_massmeanlog']    = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/gt0_massmeanlog.txt')
   dict_data['order k=%d'%(k)]['gtend_massmeanlog']  = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/gtend_massmeanlog.txt')
   
   gij_t0 = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/gij_t0.txt')
   gij_tend = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/gij_tend.txt')
   gij = np.genfromtxt(path_data+'/nbins='+str(nbins)+'/kpol='+str(k)+'/gij.txt')

   time = dict_data['order k=%d'%(k)]['time']
   if (k==0):
      gij = np.reshape(gij,(len(time),nbins))

   else:
      gij_t0 = np.reshape(gij_t0,(nbins,k+1))
      gij_tend = np.reshape(gij_tend,(nbins,k+1))
      gij = np.reshape(gij,(len(time),nbins,k+1))

   dict_data['order k=%d'%(k)]['gij_t0']   = gij_t0
   dict_data['order k=%d'%(k)]['gij_tend'] = gij_tend
   dict_data['order k=%d'%(k)]['gij']      = gij


   #compute rhodust mass fraction
   massgrid = dict_data['order k=%d'%(k)]['massgrid']

   eps_rhodust = 1e-40

   if (k==0):
      rhodust_t0 = gij_t0 * (massgrid[1:]-massgrid[:nbins])
      rhodust_tend = gij_tend * (massgrid[1:]-massgrid[:nbins])
   else:
      rhodust_t0 = gij_t0[:,0] * (massgrid[1:]-massgrid[:nbins])
      rhodust_tend = gij_tend[:,0] * (massgrid[1:]-massgrid[:nbins])

   rhodust_t0[rhodust_t0 < eps_rhodust] = eps_rhodust
   rhodust_tend[rhodust_tend < eps_rhodust] = eps_rhodust


   mass_frac_t0 = rhodust_t0/np.sum(rhodust_t0)
   mass_frac_tend = rhodust_tend/np.sum(rhodust_tend)

   dict_data['order k=%d'%(k)]['rhodust_t0']   = rhodust_t0
   dict_data['order k=%d'%(k)]['rhodust_tend'] = rhodust_tend

   dict_data['order k=%d'%(k)]['mass_fraction_t0']   = mass_frac_t0
   dict_data['order k=%d'%(k)]['mass_fraction_tend'] = mass_frac_tend

   



#define limits for plot
first_key = list(dict_data.keys())[0]
xmin = dict_data[first_key]['massgrid'][0]
xmax = dict_data[first_key]['massgrid'][-1]
tend = dict_data[first_key]['time'][-1]

ymin = 10**(-16)
ymax = 10

#mass grid for exact solution
x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000)


match kernel:
   case "kconst":
      title = "Constant kernel"
   case "kadd":
      title = "Additive kernel"
   case _:
      print("Need to add exact solution for other simple kernels")
      sys.exit()



plt.figure(1)

plt.loglog(x,exact_solution(kernel,x,0.),'--',c='C0',alpha=0.5,zorder=1)
plt.loglog(x,exact_solution(kernel,x,tend),'--',c='C0',label='Analytic',zorder=2)

for i in dict_data:
   
   kpol               = dict_data[i]['kpol']
   massgrid           = dict_data[i]['massgrid']
   massbins           = dict_data[i]['massbins']
   massmeanlog        = dict_data[i]['massmeanlog']
   gij_t0             = dict_data[i]['gij_t0']
   gij_tend           = dict_data[i]['gij_tend']
   gt0_massmeanlog    = dict_data[i]['gt0_massmeanlog']
   gtend_massmeanlog  = dict_data[i]['gtend_massmeanlog']


   if (kpol==0):
      color = 'black'
      
   else:
      color = "C"+str(kpol)

   #plot polynomials
   if (display_polynomials):
      for j in range(nbins):
         plt.plot(I(massgrid,j),recons_g(massgrid,massbins,gij_t0,kpol,j,I(massgrid,j)),c=color,alpha=0.5,zorder=3)
         plt.plot(I(massgrid,j),recons_g(massgrid,massbins,gij_tend,kpol,j,I(massgrid,j)),c=color,zorder=3)


   #plot polynomial value at geometric mean of bins
   plt.loglog(massmeanlog,gt0_massmeanlog,markeredgecolor=color,**marker_style,alpha=0.5,zorder=4)
   plt.loglog(massmeanlog,gtend_massmeanlog,markeredgecolor=color,label='order k=%d'%(kpol),**marker_style,zorder=4)



plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel(r'mass $m$')
plt.ylabel(r'mass density $g(m,\tau)$')
plt.legend(loc='lower left',ncol=1)
plt.title(title + r', $\tau=%.2f$'%(tend))
plt.tight_layout()


# to save plot
# plt.savefig('./' + kernel + '_mass_density.pdf',**savefig_options)


plt.figure(2)

plt.loglog(x,x*exact_solution(kernel,x,0.),'--',c='C0',alpha=0.5,zorder=1)
plt.loglog(x,x*exact_solution(kernel,x,tend),'--',c='C0',label='Analytic',zorder=2)



for i in dict_data:
   
   kpol               = dict_data[i]['kpol']
   massgrid           = dict_data[i]['massgrid']
   massbins           = dict_data[i]['massbins']
   massmeanlog        = dict_data[i]['massmeanlog']
   gij_t0             = dict_data[i]['gij_t0']
   gij_tend           = dict_data[i]['gij_tend']
   gt0_massmeanlog    = dict_data[i]['gt0_massmeanlog']
   gtend_massmeanlog  = dict_data[i]['gtend_massmeanlog']

   if (kpol==0):
      color = 'black'
      
   else:
      color = "C"+str(kpol)


   #plot polynomials
   if (display_polynomials):
      for j in range(nbins):
         if (kpol == 0):
            plt.plot(I(massgrid,j),massmeanlog[j]*recons_g(massgrid,massbins,gij_t0,kpol,j,I(massgrid,j)),c=color,alpha=0.5,zorder=3)
            plt.plot(I(massgrid,j),massmeanlog[j]*recons_g(massgrid,massbins,gij_tend,kpol,j,I(massgrid,j)),c=color,zorder=3)
         else:
            plt.plot(I(massgrid,j),I(massgrid,j)*recons_g(massgrid,massbins,gij_t0,kpol,j,I(massgrid,j)),c=color,alpha=0.5,zorder=3)
            plt.plot(I(massgrid,j),I(massgrid,j)*recons_g(massgrid,massbins,gij_tend,kpol,j,I(massgrid,j)),c=color,zorder=3)

   #plot polynomial value at geometric mean of bins
   plt.loglog(massmeanlog,massmeanlog*gt0_massmeanlog,markeredgecolor=color,**marker_style,alpha=0.5,zorder=4)
   plt.loglog(massmeanlog,massmeanlog*gtend_massmeanlog,markeredgecolor=color,label='order k=%d'%(kpol),**marker_style,zorder=4)



plt.xlim(xmin,xmax)
plt.ylim(1e-7,1.)
plt.xlabel(r'mass $m$')
plt.ylabel(r'distribution $m.g(m,\tau)$')

plt.legend(loc='lower left',ncol=1)
plt.title(title + r', $\tau=%.2f$'%(tend))
plt.tight_layout()


# to save plot
# plt.savefig('./' + kernel + 'mxg_distribution.pdf',**savefig_options)



#check mass conservation
#int(x*Exp[-x],{x,xmin,xmax})
def M1_exact(xmin,xmax):
   return np.exp(-xmin)*(1.+xmin)-np.exp(-xmax)*(1.+xmax)

plt.figure(3)
for i in dict_data:
   
   kpol       = dict_data[i]['kpol']
   massgrid   = dict_data[i]['massgrid']
   gij        = dict_data[i]['gij']
   time       = dict_data[i]['time']

   if (kpol==0):
      color = 'black'
      total_mass_density = np.sum(gij * (massgrid[1:]-massgrid[0:nbins]),axis=1)
   else:
      color = "C"+str(kpol)
      total_mass_density = np.sum(gij[:,:,0] * (massgrid[1:]-massgrid[0:nbins]),axis=1)

   abs_err = np.abs(total_mass_density[1:]-M1_exact(xmin,xmax))/M1_exact(xmin,xmax)
   plt.loglog(time[1:],abs_err,c=color,label='order k=%d'%(kpol))


ymin = 10**(np.floor(np.log10(np.mean(abs_err)))-2.)
ymax = 10**(np.floor(np.log10(np.mean(abs_err)))+2.)

plt.ylim(ymin,ymax)

plt.ylabel(r'Absolute error total mass density ')
plt.xlabel('time')
plt.legend(loc='lower left',ncol=1)
plt.title(title)
plt.tight_layout()


# to save plot
# plt.savefig('./' + kernel + '_mass_cons.pdf',**savefig_options)


plt.show()

