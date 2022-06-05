import argparse
import numpy as np
from latqcdtools.physics.HRG import HRG,EV_HRG
import latqcdtools.base.logger as logger

T = 156.5

parser = argparse.ArgumentParser(description='Script to carry out HRG calculations.')
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--tag", dest="particle_list", help="Name of the particle list", required=True ,type=str)
parser.add_argument('--temperature_range', dest='temperature_range',required=False, help="Perform HRG calculation in this range.",type=str)
parser.add_argument("--b", dest="b", required=True, help="Excluded volume parameter.", type=float)
#parser.add_argument('--column', dest='column',nargs='+', required=True, help="Read from these columns in HRG file.",type=int)
parser.add_argument("--obs", dest="obs", help="Observable to calculate (p, chi)", default="chi", type=str)
parser.add_argument("--bqs", dest="BQS", required=False, help="BQS mu derivative orders.", type=str)


args = parser.parse_args()    



# This is generally the QM hrg file
hadrons,M,Q,B,S,C,g=np.loadtxt(args.hadron_file.name,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")


# PDG HRG file
hadrons1,M1,Q1,B1,S1,C1,g1=np.loadtxt("../latqcdtools/physics/PDG_hadron_list_ext_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))

muB=T*np.linspace(0,3,80)

# mu derivative orders
Border = int(args.BQS[0])
Qorder = int(args.BQS[1])
Sorder = int(args.BQS[2])

b  = args.b
tag = str(args.particle_list)
# spin statistics
w  = np.array([1 if ba==0 else -1 for ba in B])
w1 = np.array([1 if ba==0 else -1 for ba in B1])


QMhrg = HRG(M,g,w,B,S,Q)
pdghrg = HRG(M1,g1,w1,B1,S1,Q1)
evhrg = EV_HRG(M,g,w,B,S,Q)
evpdghrg = EV_HRG(M1,g1,w1,B1,S1,Q1)

# initialize the class for mesons

mesons_qm=HRG(M[np.where(B==0)],g[np.where(B==0)],w[np.where(B==0)],B[np.where(B==0)],S[np.where(B==0)],Q[np.where(B==0)])
mesons_pdg=HRG(M1[np.where(B1==0)],g1[np.where(B1==0)],w1[np.where(B1==0)],B1[np.where(B1==0)],S1[np.where(B1==0)],Q1[np.where(B1==0)])

# Reminder:For Q diagonal cumulants we need pion contribution with EVHRG  
#and S diagonal cumulants we need kaon contribution with EVHRG

chi_QM = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder,mu_B=muB)
chi_pdg = pdghrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder,mu_B=muB)

if (Border==0):
    chi_ev = evhrg.gen_chi(T,b,1,Q_order=Qorder,S_order=Sorder,mu_B=muB) + \
    evhrg.gen_chi(T,b,-1,Q_order=Qorder, S_order=Sorder,mu_B=muB) + mesons_qm.gen_chi(T,Q_order=Qorder,S_order=Sorder,mu_B=muB)
    chi_ev1 = evpdghrg.gen_chi(T,b,1,Q_order=Qorder,S_order=Sorder,mu_B=muB)\
    +evpdghrg.gen_chi(T,b,-1,Q_order=Qorder,S_order=Sorder,mu_B=muB) \
    + mesons_pdg.gen_chi(T,Q_order=Qorder,S_order=Sorder,mu_B=muB)
else:
    chi_ev = evhrg.gen_chi(T,b,1, B_order=Border, Q_order=Qorder, S_order=Sorder,mu_B=muB)\
    +evhrg.gen_chi(T,b,-1, B_order=Border, Q_order=Qorder, S_order=Sorder,mu_B=muB)
    chi_ev1 = evpdghrg.gen_chi(T,b,1, B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB) \
    + evpdghrg.gen_chi(T,b,-1, B_order=Border, Q_order=Qorder, S_order=Sorder,mu_B=muB)
# 1 is for anti-particle and -1 is for anti-particle

np.savetxt("chiBQS_%s_Hrg_muB_fixefT%0.1f_f_BI_b%0.2f_%s"%(args.BQS,T,args.b,tag),
               np.c_[muB/T,chi_pdg,chi_QM,chi_ev,chi_ev1],fmt='%.1f %.8e %.8e %.8e %0.8e',
               header='muB / T    PDG-HRG         QM-HRG          EV-HRG_b%d       EV_PDGHRG_b%d' % (b, b))

