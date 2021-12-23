import math as mt
import numpy as np
from latqcdtools.tools import fm_to_GeVinv
from solve import solve


def beta_func(beta):
  b0=9/(16*np.pi**2)
  b1=1/(4*np.pi**4)
  return (((10*b0)/beta)**(-b1/(2*b0**2)))*np.exp(-beta/(20*b0))


# a/r_1 for ms/ml = 20
def a_r1(beta, error=False):
    c0=43.1     # taken from hotQCD paper
    c2=343236  # taken from hotQCD paper
    d2=5514    # taken from hotQCD paper

    res = (c0*beta_func(beta) + c2*(10/beta)*beta_func(beta)**3) / (1 + d2*(10/beta)*beta_func(beta)**2)
    
    if error:
        c0_err=0.3 # taken from hotQCD paper
        c2_err=41191 # taken from hotQCD paper
        d2_err=755 # taken from hotQCD paper

        p=1+(d2*10*beta_func(beta)**2)/beta

        a=((c0_err*beta_func(beta))/p)**2

        b=((c2_err*10*beta_func(beta)**3)/(p*beta))**2

        c=( ( -(c0*beta_func(beta) + (c2*10*beta_func(beta)**3)/beta)*(10*d2_err*beta_func(beta)**2)/beta ) / (p**2) )**2

        res_err = (mt.sqrt(a+b+c))

        return res, res_err 

    return res

# a/r_0 for ms/ml = 20
def a_r0(beta, error=False):
    r0_div_r1 = 1.5092 # taken from hotQCD paper 
    ar1, ar1_err = a_r1(beta, error=True)
    res = ar1/r0_div_r1

    if error:
        r0_div_r1_err = 0.0039 # taken from hotQCD paper
        da = ar1_err/r0_div_r1
        db = (-ar1*r0_div_r1_err)/(r0_div_r1**2)
        res_err = (mt.sqrt(da**2+db**2))
        return res, res_err 

    return res


def a_fm(beta):
    r1_phys_fm= 0.3106 # = r1_cont taken from hotQCD paper
    return r1_phys_fm*a_r1(beta)

def a_GeV(beta):
    r1_phys_fm= 0.3106 # = r1_cont taken from hotQCD paper
    return fm_to_GeVinv(r1_phys_fm*a_r1(beta))

def get_T_GeV(beta, nt):
    return 1./(nt*a_GeV(beta))

def search_beta(T_GeV, nt):
    # THE BOUDARIES ARE WRONG!!
    return solve(get_T_GeV, T_GeV, 5.9, 7.8, 10e-12, nt)