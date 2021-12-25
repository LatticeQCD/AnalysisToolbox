import numpy as np
import latqcdtools.logger as logger


def V_Teq0(r):
    """Zero temperature quark potential in [MeV], takes r in [fm]. The parameters a, b, and c come from
     a Levenberg-Marquardt fit of the data in Fig 14 of Phys. Rev. D90 (2014) 094503. These numbers can
     be obtained again by running analysistoolbox/hisq_potential/fit_hisq_pot.py."""
    #    result:  [ -91.30436191 1022.25286821  106.70659264]
    #     error:  [0.53809612 2.51598869 2.58370288]
    #  chi2/dof:  0.8083127937775374
    a=  -91.30436191
    b= 1022.25286821
    c=  106.70659264
    return a/r + b*r + c


def impdist(Ns,r2max):
    """Calculation of tree-level improved distances.

    INPUT:
         Ns--spatial extension of lattice
      r2max--maximum squared distance to improve

    OUTPUT:
      r2imp--list of improved distances"""
    r2imp=[]
    if not Ns>0:
        logger.TBError("impdist: Need Ns>0")
    kn=2.*np.pi/Ns
    pots  =[0.]*3*Ns**2
    weight=[0 ]*3*Ns**2
    cosf  =[]
    sinf  =[]
    for i in range(3*Ns**2):
        cosf.append(np.cos(i*kn))
    for i in range(Ns):
        sinf.append(np.sin(i*kn/2.))
    for x in range(int(Ns/4)+1):
        for y in range(int(Ns/4)+1):
            for z in range(int(Ns/4)+1):
                sq=x**2+y**2+z**2
                if sq>r2max:
                    continue
                pot=0.
                for k1 in range(Ns):
                    for k2 in range(Ns):
                        for k3 in range(Ns):
                            if not (k1+k2+k3)==0:
                                r1  =cosf[k1*x+k2*y+k3*z]
                                r2  =sinf[k1]**2+1./3.*sinf[k1]**4
                                r2 +=sinf[k2]**2+1./3.*sinf[k2]**4
                                r2 +=sinf[k3]**2+1./3.*sinf[k3]**4
                                pot+=r1/(4.*r2)
                pot*=1./Ns**3
                pots[sq]+=pot
                weight[sq]+=1
    for x in range(int(Ns/4)+1,int(Ns/2)+1):
        sq=x**2
        if sq>r2max:
            continue
        pot=0.
        for k1 in range(Ns):
            for k2 in range(Ns):
                for k3 in range(Ns):
                    if not (k1+k2+k3)==0:
                        r1  =cosf[k1*x]
                        r2  =sinf[k1]**2+1./3.*sinf[k1]**4
                        r2 +=sinf[k2]**2+1./3.*sinf[k2]**4
                        r2 +=sinf[k3]**2+1./3.*sinf[k3]**4
                        pot+=r1/(4.*r2)
        pot*=1./Ns**3
        pots[sq]+=pot
        weight[sq]+=1
    for i in range(1,r2max+1):
        if not weight[i]==0:
            r2imp.append(1./(4.*np.pi*(pots[i]/weight[i]+0.22578/Ns)))
    return r2imp
