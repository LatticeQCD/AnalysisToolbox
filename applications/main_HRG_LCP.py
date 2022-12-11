# 
# main_HRG_LCP.py                                                               
# 
# J. Goswami, D. Clarke
# 
# Do some HRG calculations to create a specified line of constant physics (LCP). For example you might constrain
# N_S=0 with N_Q/N_B=0.5.
#
import numpy as np
import argparse
from latqcdtools.physics.HRG import HRG, EV_HRG, LCP_init_NS0
from latqcdtools.base.utilities import getArgs, printArg
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import writeTable
from latqcdtools.math.optimize import persistentSolve
from latqcdtools.base.check import rel_check


logger.set_log_level('INFO')


parser = argparse.ArgumentParser(description='Script to determine muB, muQ and muS along the strangeness neutral trajectory',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", default="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt",
                    help="Table with hadron properties")
parser.add_argument("--b", dest="b", help="excluded volume parameter.", default=None, type=float)
parser.add_argument("--Tpc", dest="Tpc", default=None,help="determine T(muB) along the pseudo-critical line", type=float)
parser.add_argument("--T", dest="temp", default=None,help="perform all calculations at this constant T", type=float)
parser.add_argument("--r", dest="r", default=0.4, help="r=nQ/nB (RHIC is 0.4, isospin symmetric is 0.5", type=float)
parser.add_argument("--models",nargs="*",dest="models",default=['QM'],required=True,help="list of HRG models from (EV,QM) to try, default = QM",type=str)


args = getArgs(parser)


models   = args.models
b        = args.b
Tpc0     = args.Tpc
r        = args.r
temp     = args.temp
muBmax   = max(1400,18*temp)


printArg("  hadron_list:",args.hadron_file)
printArg("     b [fm^3]:",args.b)
printArg("      T [MeV]:",args.temp)
printArg(" muBmax [MeV]:",muBmax)


if "EV" in models and b is None:
    logger.TBError("Need to specify excluded volume parameter b when using EVHRG.")
if (Tpc0 is not None) and (temp is not None):
    logger.TBError("Please choose between having a fixed temperature or moving along pseudocritical line.")


muB = np.arange(10,muBmax,10)


hadrons,M,Q,B,S,C,g,w = np.loadtxt(args.hadron_file,unpack=True,dtype="U11,f8,i8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))


if Tpc0 is not None:
    kap2 = 0.012
    temp = Tpc0
    t = Tpc0*(1.0-kap2*(muB/Tpc0)**2)
else:
    t = np.repeat(temp,len(muB))


muBh  = muB/t
QMhrg = HRG(M, g, w, B, S, Q)
evhrg = EV_HRG(M, g, w, B, S, Q)


#
# strategy: given muB, solve for muQ and muS such that NS = 0 and N1/NB = r
#
def strangeness_neutral_equations(muQSh,muh,T,hrg):

    x, y   = muQSh
    maskm  = B == 0
    mesons = HRG(M[maskm],g[maskm],w[maskm],B[maskm],S[maskm],Q[maskm])

    if hrg=='QM':
        # chi_1X = N_X
        X1B = QMhrg.gen_chi(T,B_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
        X1S = QMhrg.gen_chi(T,S_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
        X1Q = QMhrg.gen_chi(T,Q_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
    else: # observables that can be improved thru EV are baryon-specific (see arXiv:2011.02812 appendix)
        X1B = evhrg.gen_chi(T,b,1,B_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y) + evhrg.gen_chi(T,b,-1, B_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
        X1Q = evhrg.gen_chi(T,b,1,Q_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y) + evhrg.gen_chi(T,b,-1, Q_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y) + mesons.gen_chi(T,Q_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
        X1S = evhrg.gen_chi(T,b,1,S_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y) + evhrg.gen_chi(T,b,-1, S_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y) + mesons.gen_chi(T,S_order=1,muB_div_T=muh,muQ_div_T=x,muS_div_T=y)
    # This is what the solver tries to make equal to zero.
    return X1S, X1Q - r*X1B


for model in models:

    print("  Solving for model",model)

    muQhi, muShi = LCP_init_NS0(muBh[0])
    muQh      = [muQhi]
    muSh      = [muShi]
    tfinal    = [t[0]]
    muBhfinal = [muBh[0]]    

    for i in range(1,len(muBh)):

        logger.details('muB/T = %8.e',muBh[i])

        muQhi, muShi = muQh[i-1], muSh[i-1]
        try:
            muQhi, muShi = persistentSolve(lambda p: strangeness_neutral_equations(p,muBh[i],t[i],model), (muQhi, muShi), tol=1e-9)
        except:
            logger.warn("No algorithm converged--giving up at muB=",muB[i])
            break

        # A quick check that the solver worked correctly.
        NB = QMhrg.gen_chi(t[i], B_order=1, muB_div_T=muBh[i], muQ_div_T=muQhi, muS_div_T=muShi)
        NQ = QMhrg.gen_chi(t[i], Q_order=1, muB_div_T=muBh[i], muQ_div_T=muQhi, muS_div_T=muShi)
        if not rel_check(NQ/NB,r,prec=1e-6):
            logger.warn('Solver discrepancy! T, NQ/NB, r = %.8e %.8e %.8e' % (t[i],NQ/NB,r))

        muQh.append(muQhi)
        muSh.append(muShi)
        tfinal.append(t[i])
        muBhfinal.append(muBh[i])

    outFileName="HRG_LCP_T%0.1f_"%temp + "r%0.1f"%r + model
    writeTable(outFileName,tfinal,muBhfinal,muQh,muSh,header=['T [MeV]','muB/T','muQ/T','muS/T'])
