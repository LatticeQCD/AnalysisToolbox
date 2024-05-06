# 
# staticPotential.py
# 
# D. Clarke
# 
# Some methods related do the static potential.
#


import numpy as np
import latqcdtools.base.logger as logger
from numba import jit


def V_Teq0(r) -> float:
    """ Zero temperature quark potential in [MeV], takes r in [fm]. The parameters a, b, and c come from
     a Levenberg-Marquardt fit of the data in Fig 14 of Phys. Rev. D90 (2014) 094503. These numbers can
     be obtained again by running analysistoolbox/hisq_potential/fit_hisq_pot.py. """
    #    result:  [ -91.30436191 1022.25286821  106.70659264]
    #     error:  [0.53809612 2.51598869 2.58370288]
    #  chi2/dof:  0.8083127937775374
    a =  -91.30436191
    b = 1022.25286821
    c =  106.70659264
    return a/r + b*r + c


def fitV_Teq0(r, a, b, c) -> float:
    """ Fit form of standard Cornell potential. Fit to be done in lattice units."""
    return a + b/r + c*r


def fitV_Teq0_oneloop(r,a,b,c,d) -> float:
    """ Including one-loop corrections to Coulomb. See Nucl. Phys. B 129 (1977) and Phys. Lett B 92 (1980).
        Fit to be done in lattice units."""
    return a + ( b + d*np.log(r) )/r + c*r


def fitV_Teq0_twoloop(r,a,b,c,d,e) -> float:
    """ Including two-loop corrections to Coulomb. See Nucl. Phys. B 501 (1997). Fit to be done in lattice units."""
    return a + ( b + d*np.log(r) + e*np.log(np.log(r)) )/r + c*r


def impdist(Ns,r2max,improvedAction=True):
    """Calculation of tree-level improved distances. Follows eq. (3) of 10.1103/PhysRevD.90.074038,

    INPUT:
           Ns--spatial extension of lattice
        r2max--maximum squared distance to improve

    OUTPUT:
        rimp--list of improved distances"""

    # This part must be placed outside the jit, since numba doesn't know what do with sys.exit.
    if not Ns > 0:
        logger.TBError("Need Ns>0")
    if r2max > (Ns/2)**2:
        logger.TBError("r2max is too large.")

    if improvedAction:
        cw = 1/3
    else:
        cw = 0

    @jit(nopython=True)
    def compiledImpDist():
        """ Ported from code by O. Kaczmarek. """
        rimp  =[]
        kn    =2.*np.pi/Ns
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
                                    r2  =sinf[k1]**2+cw*sinf[k1]**4
                                    r2 +=sinf[k2]**2+cw*sinf[k2]**4
                                    r2 +=sinf[k3]**2+cw*sinf[k3]**4
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
                            r2  =sinf[k1]**2+cw*sinf[k1]**4
                            r2 +=sinf[k2]**2+cw*sinf[k2]**4
                            r2 +=sinf[k3]**2+cw*sinf[k3]**4
                            pot+=r1/(4.*r2)
            pot*=1./Ns**3
            pots[sq]+=pot
            weight[sq]+=1
        for i in range(1,r2max+1):
            if not weight[i]==0:
                rimp.append(1./(4.*np.pi*(pots[i]/weight[i]+0.22578/Ns)))
        return rimp

    return compiledImpDist()