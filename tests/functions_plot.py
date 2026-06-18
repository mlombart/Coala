import scipy.special as special
from mpmath import *
import os
import re
import numpy as np
import sys


#legendre polynomials
def LegendreP(i,x):
   coeffs=special.legendre(i)
   res = 0
   for j in range(0,i+1):
      res += coeffs[j]*x**j

   return res

#polynomials reconstruction
def recons_g(massgrid,massbins,gij,k,j,x):
   # to map bin j onto [-1,1]
   xij = 2./(massgrid[j+1]-massgrid[j])*(x-massbins[j])
   if (k==0):
      res = np.polynomial.legendre.legval(xij, gij[j]) 
   else:
      res = np.polynomial.legendre.legval(xij, gij[j,:])   
   return res


def I(massgrid,j):
   res= np.logspace(np.log10(massgrid[j]),np.log10(massgrid[j+1]),num=200)
   return res


def get_order_kpol(path):
    try:
        list_dir = os.listdir(path=path)
        list_dir_clean = [f for f in list_dir if re.match(r'kpol.*', f)]
        list_kpol= []
        for i in list_dir_clean:
            num = re.findall('[0-9]+', i)
            list_kpol.append(int(num[0]))

    except OSError as e:
        print(f"Error occurred: {e}")
        print("Need to compute these data !")
        sys.exit()

    return sorted(list_kpol)


#solution kadd g(x,0)=x exp(-x)
def solkadd(x,tau):
    if tau==0:
        res = [x[i]*exp(-x[i]) for i in range(len(x))]
    else:
        T = 1-exp(-tau)
        res = [((1-T)*exp(-x[i]*(T+1))*besseli(1,2*x[i]*sqrt(T)))/(sqrt(T)) for i in range(np.size(x))]
    return res

#solution kconst g(x,0)=x exp(-x)
def solkconst(x,tau):
   res = 4.*x/((2+tau)**2)*np.exp(-(1-tau/(2.+tau))*x)
   return res


def exact_solution(kernel,x,tau):
    res = 0.
    match kernel:
        case "kconst":
            res = solkconst(x,tau)
        case "kadd":
            res = solkadd(x,tau)
        case _:
            print("Need to add exact solution for other simple kernels")
            sys.exit()

    return res

