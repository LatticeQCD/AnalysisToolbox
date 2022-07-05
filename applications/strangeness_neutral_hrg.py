import numpy as np
import argparse
from scipy.optimize import newton_krylov
from latqcdtools.physics.HRG import HRG,EV_HRG
import latqcdtools.base.logger as logger


parser = argparse.ArgumentParser(description='Script to determine muB, muQ and muS along the strangeness neutral trajectory',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", required=True,help="Table with hadron properties.", type=lambda f: open(f))
parser.add_argument("--tag", dest="particle_list", help="Name of the particle list", required=True ,type=str)
parser.add_argument("--b", dest="b", required=True, help="Excluded volume parameter.", type=float)
parser.add_argument("--Tpc", dest="Tpc", required=True,help="value of the pseudo-critical temperature at zero chemical potential", type=float)
parser.add_argument("--r", dest="r", required=True, help="nQ/nB = 0.4", type=float)


args, invalid_args = parser.parse_known_args()
if len(invalid_args)>0:
    logger.TBError("Received unrecognized arguments",invalid_args)

# value of the excluded parameter
b = args.b

# value of the n_Q/n_B, 0.4 = RHIC and 0.5 = Isospin symmetery
r = args.r
#This is  the hrg file
hadrons,M,Q,B,S,C,g,w=np.loadtxt(args.hadron_file.name,unpack=True,dtype="U11,f8,i8,i8,i8,i8,i8,i8",usecols=(0,1,2,3,4,5,6,7))

tag = str(args.particle_list)

maskm=B==0
mask=B==1

hadronsB=hadrons[mask]
MB=hadrons[mask]
QB=Q[mask]
BB=B[mask]
SB=S[mask]
CB=C[mask]
gB=g[mask]
wB=w[mask]

Tpc0=args.Tpc 
kap2=0.012


# This paramertization is kind of outdated however probably useful for solver to give a better starting values of muQ and muS
cf=0.0
df=1.308
ef=0.273
dS=0.214
eS=0.161
dQ=0.0211
eQ=0.106
#***************************************************************************************************#
def strangeness_neutral_equations(muQS,muB,T,hrg):
    x, y = muQS
    QMhrg  = HRG(M,g,w,B,S,Q)
    evhrg  = EV_HRG(M,g,w,B,S,Q)
    mesons = HRG(M[maskm],g[maskm],w[maskm],B[maskm],S[maskm],Q[maskm])
    if hrg=='QMhrg':
        X1B = QMhrg.gen_chi(T,B_order=1,mu_B=muB,mu_Q=x,mu_S=y)
        X1Q = QMhrg.gen_chi(T,Q_order=1,mu_B=muB,mu_Q=x,mu_S=y)
        X1S = QMhrg.gen_chi(T,S_order=1,mu_B=muB,mu_Q=x,mu_S=y)
    else:
        X1B = evhrg.gen_chi(T,b,1,B_order=1,mu_B=muB,mu_Q=x,mu_S=y)+evhrg.gen_chi(T,b,-1, B_order=1,mu_B=muB,mu_Q=x,mu_S=y)
        X1Q = evhrg.gen_chi(T,b,1,Q_order=1,mu_B=muB,mu_Q=x,mu_S=y)+evhrg.gen_chi(T,b,-1, Q_order=1,mu_B=muB,mu_Q=x,mu_S=y)+mesons.gen_chi(T,Q_order=1,mu_B=muB,mu_Q=x,mu_S=y)
        X1S = evhrg.gen_chi(T,b,1,S_order=1,mu_B=muB,mu_Q=x,mu_S=y)+evhrg.gen_chi(T,b,-1, S_order=1,mu_B=muB,mu_Q=x,mu_S=y)+mesons.gen_chi(T,S_order=1,mu_B=muB,mu_Q=x,mu_S=y)
    return X1S, X1Q - r * X1B


# Generating muB values
xs=np.arange(0.95,600,4)

#We solve for a given muB, muQ and muS, strating values of muQ and muS
fmuB =  xs
fmuQ = -dQ/(1.0 + eQ*xs) 
fmuS =  dS/(1.0 + eS*xs)




Tpc = Tpc0*(1.0-kap2*(fmuB/Tpc0)**2)


fmu_sol_qmhrg = newton_krylov(lambda p:strangeness_neutral_equations(p,fmuB,Tpc,'QMhrg'), (fmuQ, fmuS))

fmu_sol_evhrg = newton_krylov(lambda p:strangeness_neutral_equations(p,fmuB,Tpc,'EVhrg'), (fmuQ, fmuS))

fmuQ_qmhrg  = fmu_sol_qmhrg[0]
fmuS_qmhrg  = fmu_sol_qmhrg[1]
fmuQ_evhrg  = fmu_sol_evhrg[0]
fmuS_evhrg  = fmu_sol_evhrg[1]





np.savetxt("Tpc%0.1f_pseudo-muS_muB_r%0.1f%s"%(Tpc0,r,tag),np.c_[Tpc,fmuB,fmuQ_qmhrg,fmuS_qmhrg,fmuQ_evhrg,fmuS_evhrg],fmt='%0.6e %0.4f %0.6e %0.6e %0.6e %0.6e',header='Tpc muB muQ_qmhrg muS_qmhrg muQ_evhrg muS_evhrg')
