

"""
Fitting of afK and a/r1 for setting the scale of HISQ action. See the scales_hisq.py for how to calculate the temperature using the fitting parameter obtained here. The data given in #https://arxiv.org/pdf/1407.6387.pdf are corrected for LCP for afK and in small beta values for r1.

package requirements: lmfit, numpy

"""


import numpy as np
from lmfit import Model, Parameter, report_fit
from colorama import Fore, Style


C0K=7.5


x,y,y_err=np.genfromtxt('afK_LCP_corrected2021.txt',unpack=True)
x1,y1,y1_err=np.genfromtxt('a_r1_small_beta_corrected2021.txt',unpack=True,usecols=(0,1,2))

def beta_func(beta):
    """
    Beta function in lattice for finite lattice spacing
    """
    b_0=9/(16*np.pi**2)
    b_1=1/(4*np.pi**4)
    return (b_0*10/beta)**(-b_1/(2*b_0**2))*np.exp(-beta/(20*b_0))

def a_times_fk(beta,c0fk,c2fk,d2fk):
    """
    Fittnng function for afK
    """
    return (c0fk*beta_func(beta)+c2fk*10./beta*beta_func(beta)**3)/ \
            (1+d2fk*10./beta*beta_func(beta)**2)

def a_div_r1(beta,c_0,c_2,d_2):
    """
    Fitting function for a/r1
    """
    return 1/((c_0*beta_func(beta) + c_2*(10/beta)*beta_func(beta)**3) / \
    (1 + d_2*(10/beta)*beta_func(beta)**2))

gmodel = Model(a_times_fk)
gmodel_r1 = Model(a_div_r1)

N1=0
N2=15
params1 = gmodel.make_params(c2fk=1,d2fk=1)
params2 = gmodel_r1.make_params(c_0=1,c_2=1.0,d_2=1.0)
result1t = gmodel.fit(y[N1:N2], params1, c0fk=Parameter('c0fk', \
           value=C0K, vary=True),beta=x[N1:N2], weights=1/y_err[N1:N2])
result2t = gmodel_r1.fit(y1, params2, beta=x1, weights=1/y1_err)

print (Fore.RED,'afK fit results:')
print(Style.RESET_ALL)
report_fit(result1t,show_correl=False)
print ('\n')
print (Fore.RED,'r1 fit results:')
print(Style.RESET_ALL)
report_fit(result2t,show_correl=False)
print ('\n')
print (Fore.RED,'r1fK fit result:')
print(Style.RESET_ALL)
print ('r1fK=',result1t.params['c0fk'].value/result2t.params['c_0'].value)
