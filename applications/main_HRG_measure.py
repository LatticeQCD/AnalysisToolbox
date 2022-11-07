# 
# main_HRG_measure.py
# 
# J. Goswami, D. Clarke
#
# Use hadron resonance lists to measure observables. The available run modes are, with muC = 0:
#    1. Fixed muB/T, muQ = muS = 0, loop over some temperatures; or
#    2. read in a LCP file, and calculate observables at fixed T along that LCP.
#


import argparse
import numpy as np
from latqcdtools.physics.HRG import HRG
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import getArgs, printArg
from latqcdtools.base.readWrite import writeTable


parser = argparse.ArgumentParser(description='Script to carry out HRG calculations.',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", default="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt",
                    help="Table with hadron properties")
parser.add_argument("--LCP_file", dest="LCP_file", help="muB, muQ, and muS chosen to fall along some line of constant physics", type=str)
parser.add_argument('--temperature_range', dest='temperature_range',required=False, help="Perform HRG calculation in this range.",type=str)
parser.add_argument("--obs", dest="obs", help="Observable to calculate (p, chi)", default="chi", type=str)
parser.add_argument("--bqsc", dest="BQSC", required=False, help="BQSC mu derivative orders.", type=str)
parser.add_argument("--muB", dest="muBh", type=float, help="muB/T")


args              = getArgs(parser)
muBh              = args.muBh
temperature_range = args.temperature_range
LCP_file          = args.LCP_file


#
# Various checks against user error.
#
if args.obs=="chi" and args.BQSC is None:
    logger.TBError("Please specify BQSC derivative orders for chi.")
if (LCP_file is not None) and (muBh is not None):
    logger.TBError("Either specify an LCP file or muB/T, not both.")
if (LCP_file is not None) and (temperature_range is not None):
    logger.TBError("The temperatures in the LCP file will overwrite what you specify.")
if (LCP_file is None) and (muBh is None):
    logger.TBError("Either specify an LCP file or muB/T, not both.")


print("\n  observable:",args.obs)
printArg(" hadron_list:",args.hadron_file)
printArg("  BQSC deriv:",args.BQSC)
printArg("    LCP file:",LCP_file)
printArg("       muB/T:",muBh)


#
# Set up the chemical potentials based on user input.
#
if muBh is not None:
    if temperature_range is None:
        T = np.arange(4, 166, 0.5)
    else:
        t = temperature_range
        T = np.arange(float(t.split(':')[0]),float(t.split(':')[1]),float(t.split(':')[2]))
    muQh = 0.
    muSh = 0.
else:
    T, muBh, muQh, muSh = np.loadtxt(LCP_file, unpack=True, usecols=(0, 1, 2, 3))
print("     T [MeV]:",T[0],T[-1],"\n")
muCh = 0.


hadrons ,M ,Q ,B ,S ,C ,g  = np.loadtxt(args.hadron_file,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")
hadrons1,M1,Q1,B1,S1,C1,g1 = np.loadtxt("../latqcdtools/physics/HRGtables/PDG_hadron_list_ext_2020.txt",unpack=True,
                                        dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6))
w  = np.array([1 if ba==0 else -1 for ba in B])
w1 = np.array([1 if ba==0 else -1 for ba in B1])


QMhrg      = HRG(M,g,w,B,S,Q,C)
pdghrg     = HRG(M1,g1,w1,B1,S1,Q1,C1)


if args.obs == "chi":

    Border = int(args.BQSC[0])
    Qorder = int(args.BQSC[1])
    Sorder = int(args.BQSC[2])
    Corder = int(args.BQSC[3])

    chi_QM  = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, C_order=Corder,
                            muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)
    chi_pdg = pdghrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, C_order=Corder,
                             muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)

    if LCP_file is None:
        writeTable("chiBQSC%s_muBh%s.txt"%(args.BQSC,str(muBh)), T, chi_pdg, chi_QM, header=['T [MeV]','PDG-HRG','QM-HRG'])
    else:
        writeTable("chiBQSC%s_%s.txt"%(args.BQSC,LCP_file), T, chi_pdg, chi_QM, header=['T [MeV]','PDG-HRG','QM-HRG'])

elif args.obs == "p":

    p_QM = QMhrg.P_div_T4(T,muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)
    p_pdg = pdghrg.P_div_T4(T,muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)

    if LCP_file is None:
        writeTable("P_div_T4_muBh%s.txt"%str(muBh),T,p_QM,p_pdg,header=['T [MeV]','PDG-HRG','QM-HRG'])
    else:
        writeTable("P_div_T4_%s.txt"%LCP_file,T,p_QM,p_pdg,header=['T [MeV]','PDG-HRG','QM-HRG'])

else:
    logger.TBError("Observable",args.obs,"not yet implemented.")
