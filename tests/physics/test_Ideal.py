# 
# testIdeal.py                                                               
# 
# D. Clarke 
# 
# Some tests for the equation of state of an ideal gas. 
# 

import numpy as np
from latqcdtools.physics.ideal import idealGas
from latqcdtools.testing import print_results, concludeTest


def testIdeal():

    iGas = idealGas(Nf=3,Nc=3)
    
    lpass = True
    
    # There's only up to 4 powers of each chemical potential in the pressure.
    lpass *= print_results(0,iGas.gen_chi(T=1,B_order=5,Q_order=0,S_order=0,C_order=0,muB=1,muQ=0,muS=0,muC=0),
                           text='B5')
    lpass *= print_results(0,iGas.gen_chi(T=1,B_order=0,Q_order=5,S_order=0,C_order=0,muB=0,muQ=1,muS=0,muC=0),
                           text='Q5')
    lpass *= print_results(0,iGas.gen_chi(T=1,B_order=0,Q_order=0,S_order=5,C_order=0,muB=0,muQ=0,muS=1,muC=0),
                           text='S5')

    # nB at muB=0 should be zero.
    lpass *= print_results(0,iGas.gen_chi(T=1,B_order=1,Q_order=0,S_order=0,C_order=0,muB=0,muQ=0,muS=0,muC=0),
                           text='nB at muB=0')

    mu0      = 867.5309
    T0       = 436.9174
    chi11BS  = iGas.gen_chi(T=T0,B_order=1,Q_order=0,S_order=1,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    chi11BQ  = iGas.gen_chi(T=T0,B_order=1,Q_order=1,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    chi11QS  = iGas.gen_chi(T=T0,B_order=0,Q_order=1,S_order=1,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    chi2B    = iGas.gen_chi(T=T0,B_order=2,Q_order=0,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    chi2Q    = iGas.gen_chi(T=T0,B_order=0,Q_order=2,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    chi2S    = iGas.gen_chi(T=T0,B_order=0,Q_order=0,S_order=2,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    n_B      = iGas.gen_chi(T=T0,B_order=1,Q_order=0,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    n_S      = iGas.gen_chi(T=T0,B_order=0,Q_order=0,S_order=1,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    pressure = iGas.P(T=T0,muB=mu0,muQ=0,muS=0,muC=0)
    entropy  = iGas.S(T=T0,muB=mu0,muQ=0,muS=0,muC=0)
    dsdT     = iGas.d2dT2_gen_chi(T=T0,B_order=0,Q_order=0,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    s_B      = iGas.ddT_gen_chi(T=T0,B_order=1,Q_order=0,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    s_Q      = iGas.ddT_gen_chi(T=T0,B_order=0,Q_order=1,S_order=0,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    s_S      = iGas.ddT_gen_chi(T=T0,B_order=0,Q_order=0,S_order=1,C_order=0,muB=mu0,muQ=0,muS=0,muC=0)
    r        = 0.5
    a_s2     =  -chi11BS**2*chi2Q + 2*chi11BQ*chi11BS*chi11QS - chi2B*chi11QS**2 \
                -chi11BQ**2*chi2S + chi2B*chi2Q*chi2S
    b_s      = a_s2
    b_B2     = chi11QS**2 - chi2Q*chi2S
    b_Q2     = chi11BS**2 - chi2B*chi2S
    b_S2     = chi11BQ**2 - chi2B*chi2Q
    b_BQ     = 2*(chi11BQ*chi2S - chi11BS*chi11QS)
    b_BS     = 2*(chi11BS*chi2Q - chi11BQ*chi11QS)
    b_QS     = 2*(chi11QS*chi2B - chi11BQ*chi11BS)
    a_s      =  -(r**2*b_Q2 + r*b_BQ + b_B2)*n_B**2 - 0.5*(r*b_QS + b_BS)*n_B*n_S
    a_B      =  (2*b_B2 +   r*b_BQ)*n_B + 0.5*b_BS*n_S
    a_Q      =  (  b_BQ + 2*r*b_Q2)*n_B + 0.5*b_QS*n_S
    a_S      =  (  b_BS +   r*b_QS)*n_B +     b_S2*n_S
    a_B2     =  -r**2*chi2S*n_B**2 + r*chi11QS*n_B*n_S
    a_Q2     =  -chi2S*n_B**2 + chi11BS*n_B*n_S
    a_S2     = (2*r*chi11BQ - r**2*chi2B - chi2Q)*n_B**2
    a_BQ     =  2*r*chi2S*n_B**2 - (r*chi11BS + chi11QS)*n_B*n_S
    a_BS     = -2*r*(chi11QS - r*chi11BS)*n_B**2 + (  chi2Q - r*chi11BQ)*n_B*n_S
    a_QS     =  2  *(chi11QS - r*chi11BS)*n_B**2 + (r*chi2B -   chi11BQ)*n_B*n_S
    eps      = entropy*T0 - pressure + mu0*n_B
    NX       = entropy**2*a_s2 + dsdT*a_s \
                + entropy*(s_B*a_B + s_Q*a_Q + s_S*a_S) \
                + s_B**2*a_B2 + s_Q**2*a_Q2 + s_S**2*a_S2 \
                + s_B*s_Q*a_BQ + s_B*s_S*a_BS + s_Q*s_S*a_QS
    DX       = dsdT*b_s + s_B**2*b_B2 + s_Q**2*b_Q2 + s_S**2*b_S2  \
                + s_B*s_Q*b_BQ + s_B*s_S*b_BS + s_Q*s_S*b_QS
    
    # Here's a check that the isentropic speed of sound is 1/3
    cs2      = NX/( (eps + pressure)*DX )
    lpass *= print_results(1/3,cs2,text='cs2')

    # Cross check with mathematica for isospin-symmetric case 
    n_B = iGas.gen_chi(T=T0,B_order=1,Q_order=0,S_order=0,C_order=0,muB=mu0,muQ=0,muS=mu0/3,muC=0)
    n_B_math = (2/81)*mu0*(mu0**2/np.pi**2+9*T0**2)
    lpass *= print_results(n_B,n_B_math,text='nB mathematica')

    concludeTest(lpass)

 
if __name__ == '__main__':
    testIdeal()
    