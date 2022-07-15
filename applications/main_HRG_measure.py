# 
# main_HRG_measure.py
# 
# J. Goswami
#
# Use data in PDG table to carry out HRG calculations.
#


import argparse
import numpy as np
from latqcdtools.physics.HRG import HRG, EV_HRG
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import getArgs, printArg
from latqcdtools.base.readWrite import writeTable


parser = argparse.ArgumentParser(description='Script to carry out HRG calculations.',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--tag", dest="particle_list", help="Name of the particle list", default=None ,type=str)
parser.add_argument('--temperature_range', dest='temperature_range',required=False, help="Perform HRG calculation in this range.",type=str)
parser.add_argument("--obs", dest="obs", help="Observable to calculate (p, chi, cs2)", default="chi", type=str)
parser.add_argument("--bqsc", dest="BQSC", required=False, help="BQSC mu derivative orders.", type=str)
parser.add_argument("--muB", dest="muB", default=None, type=float, help="muB/T")
parser.add_argument("--models",nargs="*",dest="models",default=['QM'],required=True,help="list of HRG models from (EV,QM,PDG) to try, default = QM",type=str)
parser.add_argument("--LCP_file", dest="LCP_file", required=False, help="muB, muQ, and muS chosen to fall along some line of constant physics", type=lambda f: open(f))
parser.add_argument("--b", dest="b", help="excluded volume parameter.", default=None, type=float)


args = getArgs(parser)


models    = args.models
b         = args.b
tag       = args.particle_list
LCP_file  = args.LCP_file
muB_div_T = args.muB


#
# Various checks against user error.
#
if "EV" in models and b is None:
    logger.TBError("Need to specify excluded volume parameter b when using EVHRG.")

if (LCP_file is not None) and (muB_div_T is not None):
    logger.TBError("The LCP file already dictates the allowed possible muB.")

if args.obs=="chi" and args.BQSC is None:
    logger.TBError("Please specify BQSC derivative orders for chi.")

if args.obs=="cs2" and LCP_file is None:
    logger.TBError("c_s^2 must be calculated at N_S=0 with some fixed s/nB.")


#
# Set up the chemical potentials based on user input.
#
if LCP_file is None:
    if args.temperature_range is None:
        T = np.linspace(130, 180, 101)
    else:
        t = args.temperature_range
        T = np.arange(float(t.split(':')[0]),float(t.split(':')[1]),float(t.split(':')[2]))
    muB = muB_div_T * T
    muQ = 0.
    muS = 0.
    muC = 0.
else:
    # common parameter file for these? he somehow needs to know which model is which, the usecols vector needs
    # to be decided automatically based on that, etc
    T, muB, muQ, muS = np.loadtxt(args.LCP_file.name, unpack=True, usecols=(0, 1, 2, 3))
    muC = 0.


print("\n  observable:",args.obs)
printArg("  BQSC deriv:",args.BQSC)
printArg("       muB/T:",muB_div_T)
printArg("    LCP file:",LCP_file.name)
printArg("         tag:",tag)
print("     T [MeV]:",T[0],T[-1],"\n")


hadrons,M,Q,B,S,C,g = np.loadtxt(args.hadron_file.name,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")
hadrons1,M1,Q1,B1,S1,C1,g1=np.loadtxt("../latqcdtools/physics/HRGtables/PDG_hadron_list_ext_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))
w  = np.array([1 if ba==0 else -1 for ba in B])
w1 = np.array([1 if ba==0 else -1 for ba in B1])



QMhrg      = HRG(M,g,w,B,S,Q,C)
pdghrg     = HRG(M1,g1,w1,B1,S1,Q1,C1)



if args.obs == "chi":

    # mu derivative orders
    Border = int(args.BQSC[0])
    Qorder = int(args.BQSC[1])
    Sorder = int(args.BQSC[2])
    Corder = int(args.BQSC[3])

    chi_QM     = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, C_order=Corder,
                               mu_B=muB, mu_Q=muQ, mu_S=muS, mu_C=muC)
    chi_pdg    = pdghrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, C_order=Corder,
                                mu_B=muB, mu_Q=muQ, mu_S=muS, mu_C=muC)

    if tag is not None:
        writeTable("chiBQSC_%s_%s.txt"%(args.BQSC,tag), T, muB_div_T, chi_pdg, chi_QM, header='T    PDG-HRG         QM-HRG  ' )
    else:
        writeTable("chiBQSC_%s_.txt"%args.BQSC, T, muB_div_T, chi_pdg, chi_QM, header='T    PDG-HRG         QM-HRG  ' )

elif args.obs == "p":

    p_QM = QMhrg.P_div_T4(T,mu_B=muB, mu_Q=muQ, mu_S=muS, mu_C=muC)
    p_pdg = pdghrg.P_div_T4(T,mu_B=muB, mu_Q=muQ, mu_S=muS, mu_C=muC)

    if tag is not None:
        writeTable("P_div_T4_%s.txt"%tag, T, p_QM, p_pdg, header='T    PDG-HRG         QM-HRG  ')
    else:
        writeTable("P_div_T4.txt",T,p_QM,p_pdg,header='T    PDG-HRG         QM-HRG  ')

elif args.obs == "cs2":

    def cs2(self, T, mu_B=0., mu_S=0., mu_Q=0., mu_C=0.):
        return (4 * self.P_div_T4(T, mu_B, mu_S, mu_Q, mu_C) + T * self.ddT_P_div_T4(T, mu_B, mu_S, mu_Q, mu_C)
                ) / (4 * self.E_div_T4(T, mu_B, mu_S, mu_Q, mu_C) + T * self.ddT_E_div_T4(T, mu_B, mu_S, mu_Q, mu_C))


