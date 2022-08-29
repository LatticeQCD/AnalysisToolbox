# 
# main_HRG_measure.py
# 
# J. Goswami
#
# Use hardon resonance lists to measure observables. You can specify some line of constant physics as input, or
# you can specify the values of the control parameters by hand.
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
parser.add_argument('--temperature_range', dest='temperature_range',required=False, help="Perform HRG calculation in this range.",type=str)
parser.add_argument("--obs", dest="obs", help="Observable to calculate (p, chi, cs2)", default="chi", type=str)
parser.add_argument("--bqsc", dest="BQSC", required=False, help="BQSC mu derivative orders.", type=str)
parser.add_argument("--muB", dest="muBh", required=True, type=float, help="muB/T")


args = getArgs(parser)
muBh = args.muBh


#
# Various checks against user error.
#
if args.obs=="chi" and args.BQSC is None:
    logger.TBError("Please specify BQSC derivative orders for chi.")


#
# Set up the chemical potentials based on user input.
#
if args.temperature_range is None:
    T = np.linspace(3, 160, 158)
else:
    t = args.temperature_range
    T = np.arange(float(t.split(':')[0]),float(t.split(':')[1]),float(t.split(':')[2]))
muQh = 0.
muSh = 0.
muCh = 0.


print("\n  observable:",args.obs)
printArg(" hadron_list:",args.hadron_file)
printArg("  BQSC deriv:",args.BQSC)
printArg("       muB/T:",muBh)
print("     T [MeV]:",T[0],T[-1],"\n")


hadrons,M,Q,B,S,C,g = np.loadtxt(args.hadron_file,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")
hadrons1,M1,Q1,B1,S1,C1,g1=np.loadtxt("../latqcdtools/physics/HRGtables/PDG_hadron_list_ext_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))
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

    writeTable("chiBQSC_%s.txt"%args.BQSC, T, muBh, chi_pdg, chi_QM, header='T    PDG-HRG         QM-HRG  ' )

elif args.obs == "p":

    p_QM = QMhrg.P_div_T4(T,muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)
    p_pdg = pdghrg.P_div_T4(T,muB_div_T=muBh, muQ_div_T=muQh, muS_div_T=muSh, muC_div_T=muCh)

    writeTable("P_div_T4.txt",T,p_QM,p_pdg,header='T    PDG-HRG         QM-HRG  ')

else:
    logger.TBError("Observable",args.obs,"not yet implemented.")
