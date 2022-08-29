# calculate susceptibilities of arbitrary order also on PC line, uses main_HRG as input, since the krylov solver
# part is the slow part we solve only once

import numpy as np
import argparse
from latqcdtools.physics.HRG import HRG
from latqcdtools.base.utilities import getArgs


parser = argparse.ArgumentParser(description='Script to calculate chiBQS along the pseudo-critical line',allow_abbrev=False)
parser.add_argument("--hadron_file", dest="hadron_file", default="../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt",
                    help="Table with hadron properties")
parser.add_argument("--fixedmuBNszerofile", dest="fixedmuB", required=True,help="values of muB, muQ and muS for fixed muB / T", type=lambda f: open(f))
parser.add_argument("--tag", dest="particle_list", help="Name of the particle list", required=True ,type=str)
parser.add_argument("--bqs", dest="BQS", required=False, help="BQS mu derivative orders.", type=str)
parser.add_argument("--obs", dest="obs", help="Observable to calculate (pressure, chi, energy , specificheat, cs2)", default="chi", type=str)
parser.add_argument("--r", dest="r", required=True, help="nQ/nB = 0.4", type=float)


args = getArgs(parser)
r = args.r

hadrons,M,Q,B,S,C,g = np.loadtxt(args.hadron_file,unpack=True,usecols=(0,1,2,3,4,5,6),dtype="U11,f8,i8,i8,i8,i8,i8")

tag = str(args.particle_list)

# spin statistics
w  = np.array([1 if ba==0 else -1 for ba in B])

#muB muQ muS file
tpc, muBh, muQh, muSh = np.loadtxt(args.fixedmuB.name,unpack=True,usecols=(0,1,2,3))

T=tpc


QMhrg = HRG(M,g,w,B,S,Q)

if args.obs == "chi":
    # mu derivative orders
    Border = int(args.BQS[0])
    Qorder = int(args.BQS[1])
    Sorder = int(args.BQS[2])
    obs = QMhrg.gen_chi(T, B_order=Border, Q_order=Qorder, S_order=Sorder, muB_div_T = muB, muS_div_T = muS, muQ_div_T = muQ)
    np.savetxt("muB_divT%0.1f_%sBQS%s_Hrg_BI_%sr%0.1f" % (muB[0]/T[0],args.obs,args.BQS,tag,r), np.c_[T,muB,obs], fmt='%.4f %0.4e %0.6e', header='T muB HRG')
elif args.obs == "pressure":
    obs = QMhrg.P_div_T4(T, muB_div_T=muB, muS_div_T = muS, muQ_div_T = muQ)
    np.savetxt("muB_divT%0.1f_%s_Hrg_BI_%sr%0.1f" % (muB[0]/T[0],args.obs,tag,r), np.c_[T,muB,obs], fmt='%.4f %0.4e %0.6e', header='T muB HRG')
elif args.obs == "energy":
    obs = QMhrg.E_div_T4(T, muB_div_T=muB, muS_div_T = muS, muQ_div_T = muQ)
    np.savetxt("muB_divT%0.1f_%s_Hrg_BI_%sr%0.1f" % (muB[0]/T[0],args.obs,tag,r), np.c_[T,muB,obs], fmt='%.4f %0.4e %0.6e', header='T muB HRG')
elif args.obs == "specificheat":
    obs = QMhrg.CV_div_T3(T, muB_div_T=muB, muS_div_T = muS, muQ_div_T = muQ)
    np.savetxt("muB_divT%0.1f_%s_Hrg_BI_%sr%0.1f" % (muB[0]/T[0],args.obs,tag,r), np.c_[T,muB,obs], fmt='%.4f %0.4e %0.6e', header='T muB HRG')
elif args.obs == "cs2":
    num = 4*QMhrg.P_div_T4(T, muB_div_T=muB,muS_div_T = muS, muQ_div_T = muQ) + T * QMhrg.ddT_P_div_T4(T, muB_div_T=muB,muS_div_T = muS, muQ_div_T = muQ)
    den = 4*QMhrg.E_div_T4(T, muB_div_T=muB,muS_div_T = muS, muQ_div_T = muQ) + T * QMhrg.ddT_E_div_T4(T, muB_div_T=muB,muS_div_T = muS, muQ_div_T = muQ)
    np.savetxt("muB_divT%0.1f_%s_Hrg_BI_%sr%0.1f" % (muB[0]/T[0],args.obs,tag,r), np.c_[T,muB,num/den], fmt='%.4f %0.4e %0.6e', header='T muB HRG')

