# 
# testHRG.py                                                               
# 
# D. Clarke 
# 
# A test for some of the HRG methods. 
#

import numpy as np
from latqcdtools.physics.HRG import HRG,EV_HRG
from latqcdtools.base.check import print_results

EPSILON = 1e-6

T = np.linspace(130, 179.5, 100)

# QM and PDG HRG files
hadrons,M,Q,B,S,C,g=np.loadtxt("../../latqcdtools/physics/QM_hadron_list_ext_strange_2020.txt",unpack=True,
                               usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")
hadrons1,M1,Q1,B1,S1,C1,g1=np.loadtxt("../../latqcdtools/physics/PDG_hadron_list_ext_2020.txt",unpack=True,
                                      dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))

# spin statistics
w  = np.array([1 if ba==0 else -1 for ba in B])
w1 = np.array([1 if ba==0 else -1 for ba in B1])

QMhrg      = HRG(M,g,w,B,S,Q)
pdghrg     = HRG(M1,g1,w1,B1,S1,Q1)
evhrg      = EV_HRG(M,g,w,B,S,Q)
evpdghrg   = EV_HRG(M1,g1,w1,B1,S1,Q1)
mesons_qm  = HRG(M[np.where(B==0)]  ,g[np.where(B==0)]  ,w[np.where(B==0)]  ,B[np.where(B==0)]  ,S[np.where(B==0)]  ,Q[np.where(B==0)])
mesons_pdg = HRG(M1[np.where(B1==0)],g1[np.where(B1==0)],w1[np.where(B1==0)],B1[np.where(B1==0)],S1[np.where(B1==0)],Q1[np.where(B1==0)])

Border    = 2
Qorder    = 0
Sorder    = 0
b         = 1
muB_div_T = 1
muB       = muB_div_T * T

chi_QM  = QMhrg.gen_chi(T , B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB)
chi_pdg = pdghrg.gen_chi(T, B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB)
if Border == 0:
    chi_ev  = evhrg.gen_chi(T, b, 1, Q_order=Qorder, S_order=Sorder, mu_B=muB) \
              + evhrg.gen_chi(T, b, -1, Q_order=Qorder, S_order=Sorder, mu_B=muB) \
              + mesons_qm.gen_chi(T, Q_order=Qorder, S_order=Sorder)
    chi_ev1 = evpdghrg.gen_chi(T, b, 1, Q_order=Qorder, S_order=Sorder, mu_B=muB) \
              + evpdghrg.gen_chi(T, b, -1, Q_order=Qorder, S_order=Sorder, mu_B=muB) \
              + mesons_pdg.gen_chi(T, Q_order=Qorder, S_order=Sorder, mu_B=muB)
else:
    chi_ev = evhrg.gen_chi(T, b, 1, B_order=Border, Q_order=Qorder, S_order=Sorder) \
             + evhrg.gen_chi(T, b, -1, B_order=Border, Q_order=Qorder, S_order=Sorder)
    chi_ev1 = evpdghrg.gen_chi(T, b, 1, B_order=Border, Q_order=Qorder, S_order=Sorder) \
              + evpdghrg.gen_chi(T, b, -1, B_order=Border, Q_order=Qorder, S_order=Sorder)

#           np.c_[T, chi_pdg, chi_QM, chi_ev, chi_ev1]

refT, refPDG, refQM, refEV, refEV1 = np.loadtxt("chiBQS_200_muB1.00_b1.00_QMHRG2020_BI.control",unpack=True)

print_results(T      , refT  , prec=EPSILON, text="T check")
print_results(chi_pdg, refPDG, prec=EPSILON, text="chiB2 PDG check")
print_results(chi_QM , refQM , prec=EPSILON, text="chiB2 QM check")
print_results(chi_ev , refEV , prec=EPSILON, text="chiB2 EV check")
print_results(chi_ev1, refEV1, prec=EPSILON, text="chiB2 EV1 check")
