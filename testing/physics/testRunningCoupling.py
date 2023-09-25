# 
# testRunningCoupling.py                                                               
# 
# D. Clarke 
# 
# Test some general methods for QCD running coupling. 
# 


import numpy as np
from latqcdtools.physics.runningCoupling import CA, CF, ZETA_3, b2_dimreg_MSbar, b3_dimreg_MSbar
from latqcdtools.math.math import print_results
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testRunningCoupling():
    
    Nf=3
    Nc=3

    logger.info('SU(3) tests:')

    print_results(CA(Nc),Nc,text='C_A from b0')
    print_results(34*CA(Nc)**2/3,102,text='C_A from b1')
    print_results(2*CF(Nc)+10*CA(Nc)/3,38/3,text='C_F from b1')

    TEST=2857./2.-5033./18.*Nf+325./54.*Nf**2
    print_results(TEST,b2_dimreg_MSbar(Nf,Nc),text='b2 MS-bar dim. reg.')

    TEST=149753./6.+3564*ZETA_3-(1078361./162.+6508./27.*ZETA_3)*Nf+(50065./162.+6472./81.*ZETA_3)*Nf**2+1093./729.*Nf**3
    print_results(TEST,b3_dimreg_MSbar(Nf,Nc),text='b3 MS-bar dim. reg.')


if __name__ == '__main__':
    testRunningCoupling()