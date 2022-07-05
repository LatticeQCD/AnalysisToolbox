import numpy as np
import argparse
from latqcdtools.physics.HRG import HRG,EV_HRG
import latqcdtools.base.logger as logger


parser = argparse.ArgumentParser(description='Script to calculate chiBQS along the pseudo-critical line',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--mubmus_file", dest="mubmus_file", required=True,help="values of muB, muQ and muS for the strangeness neutral case in the pseudo-critical line", type=lambda f: open(f))
parser.add_argument("--tag", dest="particle_list", help="Name of the particle list", required=True ,type=str)
parser.add_argument("--bqs", dest="BQS", required=False, help="BQS mu derivative orders.", type=str)
parser.add_argument("--b", dest="b", required=True, help="Excluded volume parameter.", type=float)
parser.add_argument("--r", dest="r", required=True, help="nQ/nB = 0.4", type=float)


args, invalid_args = parser.parse_known_args()
if len(invalid_args)>0:
    logger.TBError("Received unrecognized arguments",invalid_args)


# This is generally the QMHRG file
hadrons,M,Q,B,S,C,g=np.loadtxt(args.hadron_file.name,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")


tag = str(args.particle_list)


# spin statistics
w  = np.array([1 if ba==0 else -1 for ba in B])

#muB muQ muS file
tpc, muB, muQ, muS, muQ_ev, muS_ev = np.loadtxt(args.mubmus_file.name,unpack=True,usecols=(0,1,2,3,4,5))

T=tpc


# mu derivative orders
Border = int(args.BQS[0])
Qorder = int(args.BQS[1])
Sorder = int(args.BQS[2])

b   = args.b

# initialize the class
QMhrg = HRG(M,g,w,B,S,Q)
evhrg = EV_HRG(M,g,w,B,S,Q)

# initialize the class for mesons

mesons_qm=HRG(M[np.where(B==0)],g[np.where(B==0)],w[np.where(B==0)],B[np.where(B==0)],S[np.where(B==0)],Q[np.where(B==0)])

# Reminder:For Q diagonal cumulants we need pion contribution with EVHRG  
#and S diagonal cumulants we need kaon contribution with EVHRG


chi_QM = QMhrg.gen_chi(T,B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B = muB, mu_S = muS, mu_Q = muQ)

if Border==0:
    chi_ev = evhrg.gen_chi(T,b,1,Q_order=Qorder,S_order=Sorder, mu_B = muB, mu_S = muS_ev, mu_Q = muQ_ev ) + evhrg.gen_chi(T,b,-1,Q_order=Qorder, S_order=Sorder , mu_B = muB, mu_S = muS_ev, mu_Q = muQ_ev) + mesons_qm.gen_chi(T, Q_order=Qorder, S_order=Sorder , mu_B = muB, mu_S = muS_ev, mu_Q = muQ_ev)
else:
    chi_ev = evhrg.gen_chi(T,b,1, B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B = muB, mu_S = muS_ev, mu_Q = muQ_ev ) + evhrg.gen_chi(T,b,-1, B_order=Border, Q_order=Qorder, S_order=Sorder, mu_B = muB, mu_S = muS_ev, mu_Q = muQ_ev)

np.savetxt("T%0.1f_pseudo-chiBQS_%s_Hrg_BI_b%0.1f%s" % (T[0],args.BQS,args.b,tag), np.c_[T,muB,chi_QM,chi_ev], fmt='%.4f %0.4e %.8e %.8e', header='T muB HRG EV-HRG_b%0.1f' % b)

