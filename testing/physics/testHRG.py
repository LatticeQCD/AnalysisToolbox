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
from latqcdtools.base.plotting import plot_lines,plot_file,set_params,latexify
import matplotlib.pyplot as plt

EPSILON = 1e-6

T = np.linspace(130, 179.5, 100)

# QM and PDG HRG files
hadrons,M,Q,B,S,C,g = np.loadtxt("../../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt",unpack=True,
                                 usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")
hadrons1,M1,Q1,B1,S1,C1,g1 = np.loadtxt("../../latqcdtools/physics/HRGtables/PDG_hadron_list_ext_2020.txt",unpack=True,
                                        dtype="U11,f8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))

# Spin statistics. w= fermi(-1)/bose(1) statistics. (If baryon number is 1, you have three spin-1/2 constituents, which
# allows only half-integer baryon spins.
w  = np.array([1 if ba==0 else -1 for ba in B])
w1 = np.array([1 if ba==0 else -1 for ba in B1])

QMhrg      = HRG(M,g,w,B,S,Q)
pdghrg     = HRG(M1,g1,w1,B1,S1,Q1)
evhrg      = EV_HRG(M,g,w,B,S,Q)
evpdghrg   = EV_HRG(M1,g1,w1,B1,S1,Q1)
mesons_qm  = HRG(M[np.where(B==0)]  ,g[np.where(B==0)]  ,w[np.where(B==0)]  ,B[np.where(B==0)]  ,S[np.where(B==0)]  ,Q[np.where(B==0)])
mesons_pdg = HRG(M1[np.where(B1==0)],g1[np.where(B1==0)],w1[np.where(B1==0)],B1[np.where(B1==0)],S1[np.where(B1==0)],Q1[np.where(B1==0)])


#
# Test: Calculate chi^200_BQS with b=1 and mu/T=1. Compare against trusted control result.
#
Border    = 2
Qorder    = 0
Sorder    = 0
b         = 1
muB_div_T = 1
muB       = muB_div_T * T

chi_QM  = QMhrg.gen_chi(T , B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB)
chi_pdg = pdghrg.gen_chi(T, B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B=muB)
chi_ev  = evhrg.gen_chi(T, b, 1, B_order=Border, Q_order=Qorder, S_order=Sorder) \
          + evhrg.gen_chi(T, b, -1, B_order=Border, Q_order=Qorder, S_order=Sorder)
chi_ev1 = evpdghrg.gen_chi(T, b, 1, B_order=Border, Q_order=Qorder, S_order=Sorder) \
          + evpdghrg.gen_chi(T, b, -1, B_order=Border, Q_order=Qorder, S_order=Sorder)

refT, refPDG, refQM, refEV, refEV1 = np.loadtxt("HRGcontrol/chiBQS_200_muB1.00_b1.00_QMHRG2020_BI.control",unpack=True)

print_results(T      , refT  , prec=EPSILON, text="T check")
print_results(chi_pdg, refPDG, prec=EPSILON, text="chiB2 PDG check")
print_results(chi_QM , refQM , prec=EPSILON, text="chiB2 QM check")
print_results(chi_ev , refEV , prec=EPSILON, text="chiB2 EV check")
print_results(chi_ev1, refEV1, prec=EPSILON, text="chiB2 EV1 check")


#
# Test: Calculate the pressure and compare it with PHYSICAL REVIEW D 90, 094503 (2014)
#
refT, ref3p_div_T4 = np.loadtxt("HRGcontrol/2014_3P_div_T4.d",unpack=True)
test3p_div_T4      = 3*pdghrg.pressure(refT,0,0,0)
print_results(ref3p_div_T4, test3p_div_T4, prec=2e-2, text="2014 HotQCD 3p/T^4 check")

#latexify()
#plot_lines(T,3*pdghrg.pressure(T,0,0,0),xmax=175,label="$\\mu_B/T=1$",marker=None)
#plot_file("HRGcontrol/2014_3P_div_T4.d",style="lines",marker=None,label="HotQCD 2014")
#set_params(xlabel="$T$",ylabel="$3P/T^4$")
#plt.show()