# 
# main_evhrg.py                                                               
# 
# J. Goswami
#
# Use data in PDG table to carry out HRG calculations.
#


import argparse
import numpy as np
from latqcdtools.physics.HRG import HRG_charm
import latqcdtools.base.logger as logger


parser = argparse.ArgumentParser(description='Script to carry out HRG calculations.')
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--tag", dest="particle_list", help="Name of the particle list", required=True ,type=str)
parser.add_argument('--temperature_range', dest='temperature_range',required=False, help="Perform HRG calculation in this range.",type=str)
parser.add_argument("--obs", dest="obs", help="Observable to calculate (p, chi)", default="chi", type=str)
parser.add_argument("--bqsc", dest="BQSC", required=False, help="BQSC mu derivative orders.", type=str)
parser.add_argument("--muB", dest="muB", default=0.0, type=float, help="muB/T")


args      = parser.parse_args()
muB_div_T = float(args.muB)


if args.obs=="chi" and args.BQSC is None:
  logger.TBError("Please specify BQSC derivative orders for chi.")


if args.temperature_range is None:
    T = np.linspace(130, 180, 101)
else:
    t = args.temperature_range
    T = np.arange(float(t.split(':')[0]),float(t.split(':')[1]),float(t.split(':')[2]))


print("\n  observable:",args.obs)
if args.BQSC is not None:
    print("   BQSC deriv:",args.BQSC)
if muB_div_T is not None:
    print("       muB/T:",muB_div_T)
print("     T [MeV]:",T[0],T[-1],"\n")


# hadron  = Name of the hadron [All particles and anti-particles]
# M       = Mass of the hadron
# Q       = charge of the hadron
# B       = Baryon number of the hadron [B=1,-1,0,0 for baryon,anti-baryon,meson,anti-mesons]
# S       = Strangeness number of the hadron
# C       = Charm number of the hadron
# g       = degenracy of the hadron state


# This is generally the QM hrg file
hadrons,M,Q,B,S,C,g=np.loadtxt(args.hadron_file.name,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")


tag = str(args.particle_list)


# PDG HRG file
hadrons1,M1,Q1,B1,S1,C1,g1=np.loadtxt("../latqcdtools/physics/PDG_hadron_list_ext_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))


# spin statistics
w  = np.array([1 if ba==0 else -1 for ba in B])
w1 = np.array([1 if ba==0 else -1 for ba in B1])


muB = muB_div_T * T


QMhrg      = HRG_charm(M,g,w,B,S,Q,C)
pdghrg     = HRG_charm(M1,g1,w1,B1,S1,Q1,C1)


if args.obs == "chi":

    # mu derivative orders
    Border = int(args.BQSC[0])
    Qorder = int(args.BQSC[1])
    Sorder = int(args.BQSC[2])
    Corder = int(args.BQSC[3])

    chi_QM     = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB)
    chi_pdg    = pdghrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB)
    # Save the output
    np.savetxt("chiBQSC_%s_muB%0.2f_%s.txt"%(args.BQSC,muB_div_T,tag),
               np.c_[T,chi_pdg,chi_QM],fmt='%.1f %.8e %.8e ',
               header='T    PDG-HRG         QM-HRG  ' )

elif args.obs == "p":
    p_QM = QMhrg.pressure(T, mu_B=muB)
    p_pdg = pdghrg.pressure(T, mu_B=args.muB)
    np.savetxt("pressure_muB%0.2f_%s.txt"%(muB_div_T,tag), np.c_[T,p_pdg,p_QM],fmt='%.1f %.8e %.8e',
               header='T    PDG-HRG         QM-HRG  ')

