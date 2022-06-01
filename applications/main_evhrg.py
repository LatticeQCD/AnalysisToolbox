# 
# main_evhrg.py                                                               
# 
# J. Goswami
#
# Use data in PDG table to carry out HRG calculations.
#


import argparse
import numpy as np
from latqcdtools.physics.HRG import HRG,EV_HRG


parser = argparse.ArgumentParser(description='Script to carry out HRG calculations.')
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--bqs", dest="BQS", required=True,help="BQS mu derivative orders.", type=str)
parser.add_argument("--b", dest="b", required=True, help="Excluded volume parameter.", type=float)
parser.add_argument('--column', dest='test4',nargs='+', required=True, help="test4",type=int)
parser.add_argument('--temperature_range', dest='temperature_range',required=False, help="Perform HRG calculation in this range.",type=str)

args = parser.parse_args()
col_tuple = tuple(args.test4)

if args.temperature_range is None:
    T=np.linspace(130, 180, 101)
else:
    t= args.temperature_range
    T=np.arange(float(t.split(':')[0]),float(t.split(':')[1]),float(t.split(':')[2]))


# hadron  = Name of the hadron [All particles and anti-particles]
# M       = Mass of the hadron
# Q       = charge of the hadron
# B       = Baryon number of the hadron [B=1,-1,0,0 for baryon,anti-baryon,meson,anti-mesons]
# S       = Strangeness number of the hadron
# C       = Charm number of the hadron
# g       = degenracy of the hadron state


#This is generally the QM hrg file
hadrons,M,Q,B,S,C,g=np.loadtxt(args.hadron_file.name,unpack=True,usecols=col_tuple,dtype="U11,f8,i8,i8,i8,i8,i8")
tag='QM'

# PDG HRG file
hadrons1,M1,Q1,B1,S1,C1,g1=np.loadtxt("../latqcdtools/physics/PDG_hadron_list_ext_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))


#spin statistics
w=np.array([1 if ba==0 else -1 for ba in B])
w1=np.array([1 if ba==0 else -1 for ba in B1])

# mu derivative orders
Border = int(args.BQS[0])
Qorder = int(args.BQS[1])
Sorder = int(args.BQS[2])

b = args.b

# initialize the class
QMhrg    = HRG(M,g,w,B,S,Q)
pdghrg   = HRG(M1,g1,w1,B1,S1,Q1)
evhrg    = EV_HRG(M,g,w,B,S,Q)
evpdghrg = EV_HRG(M1,g1,w1,B1,S1,Q1)

# initialize the class for mesons

mesons_qm=HRG(M[np.where(B==0)],g[np.where(B==0)],w[np.where(B==0)],B[np.where(B==0)],S[np.where(B==0)],Q[np.where(B==0)])
mesons_pdg=HRG(M1[np.where(B1==0)],g1[np.where(B1==0)],w1[np.where(B1==0)],B1[np.where(B1==0)],S1[np.where(B1==0)],Q1[np.where(B1==0)])

# Reminder:For Q diagonal cumulants we need pion contribution with EVHRG
#and S diagonal cumulants we need kaon contribution with EVHRG

chi_QM = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder)
chi_pdg = pdghrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder)

if Border==0:
    chi_ev = evhrg.gen_chi(T,b,1,Q_order=Qorder,S_order=Sorder)+evhrg.gen_chi(T,b,-1,Q_order=Qorder, S_order=Sorder) + mesons_qm.gen_chi(T,Q_order=Qorder,S_order=Sorder)
    chi_ev1 = evpdghrg.gen_chi(T,b,1,Q_order=Qorder,S_order=Sorder)+evpdghrg.gen_chi(T,b,-1,Q_order=Qorder,S_order=Sorder) + mesons_pdg.gen_chi(T,Q_order=Qorder,S_order=Sorder)
else:
    chi_ev = evhrg.gen_chi(T,b,1, B_order=Border, Q_order=Qorder, S_order=Sorder)+evhrg.gen_chi(T,b,-1, B_order=Border, Q_order=Qorder, S_order=Sorder)
    chi_ev1 = evpdghrg.gen_chi(T,b,1, B_order=Border, Q_order=Qorder, S_order=Sorder) + evpdghrg.gen_chi(T,b,-1, B_order=Border, Q_order=Qorder, S_order=Sorder)

# 1 is for anti-particle and -1 is for anti-particle

#Save the output
filename='/'.join(args.hadron_file.name.split('/')[1:])
np.savetxt("chiBQS_%s_Hrg_BI_b%0.2f_%s"%(args.BQS,args.b,tag),
           np.c_[T,chi_pdg,chi_QM,chi_ev,chi_ev1],fmt='%.1f %.8e %.8e %.8e %0.8e',
           header='T PDG-HRG QM-HRG EV-HRG_b%d EV_PDGHRG_b%d'%(b,b))

print ("Done. Check the data files chiBQS*")

