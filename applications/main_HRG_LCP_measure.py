# 
# main_HRG_cs2.py
# 
# D. Clarke
#
# Do an HRG calculation of basic observables along a line of constant entropy.
#


import argparse
import numpy as np
from latqcdtools.physics.HRG import HRG, LCP_init_NS0
from latqcdtools.base.utilities import getArgs, find_nearest_idx
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, newton_krylov
from latqcdtools.base.check import rel_check


parser = argparse.ArgumentParser(description='Script to calculate cs2 in HRG.',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--LCP_file", dest="LCP_file", required=True, help="muB, muQ, and muS chosen to fall along some line of constant physics", type=lambda f: open(f))
parser.add_argument("--snB", dest="snB", required=True, help="fix s/nB to this value", type=float)
parser.add_argument("--r", dest="r", default=0.4, help="r=nQ/nB (RHIC is 0.4, isospin symmetric is 0.5", type=float)


args = getArgs(parser)

showPlots  = False
target_snB = args.snB
LCP_file   = args.LCP_file
r          = args.r


T, muBh, muQh, muSh = np.loadtxt(args.LCP_file.name, unpack=True, usecols=(0, 1, 2, 3))
muCh = 0.


hadrons,M,Q,B,S,C,g = np.loadtxt(args.hadron_file.name,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")
w  = np.array([1 if ba==0 else -1 for ba in B])

QMhrg      = HRG(M,g,w,B,S,Q,C)

s  = QMhrg.S_div_T3(T,muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)
nB = QMhrg.gen_chi(T, B_order=1, Q_order=0, S_order=0, C_order=0, muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)

x=muBh
y=s/nB

guessIndex=find_nearest_idx(y,target_snB)

yinterp = interp1d(x, y, kind='cubic')
def func(z):
    return yinterp(z) - target_snB
target_muBh = fsolve(func,x[guessIndex])[0]


def strangeness_neutral_equations(muQSh,muh,t):
    x, y   = muQSh
    X1B = QMhrg.gen_chi(t,B_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
    X1S = QMhrg.gen_chi(t,S_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
    X1Q = QMhrg.gen_chi(t,Q_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
    return X1S, X1Q - r*X1B


muQh, muSh = LCP_init_NS0(muBh)
solution = newton_krylov(lambda p: strangeness_neutral_equations(p,target_muBh,T[0]), (muQh, muSh))
muQh = solution[0]
muSh = solution[1]


s  = QMhrg.S_div_T3(T[0],muB_div_T=target_muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)
nB = QMhrg.gen_chi(T[0], B_order=1, Q_order=0, S_order=0, C_order=0, muB_div_T=target_muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)


if not rel_check(target_snB,s/nB,1e-4):
    print("#Strayed from line of constant physics by",100*abs(target_snB-s/nB)/target_snB,"%")

pT4 = QMhrg.P_div_T4(T[0],muB_div_T=target_muBh,muS_div_T=muSh,muQ_div_T=muQh,muC_div_T=muCh)
eT4 = QMhrg.E_div_T4(T[0],muB_div_T=target_muBh,muS_div_T=muSh,muQ_div_T=muQh,muC_div_T=muCh)


print('  %.8e  %.8e  %.8e  %.8e'%(T[0],target_muBh,pT4,eT4))
