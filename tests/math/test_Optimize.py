# 
# test_Optimize.py                                                               
# 
# D. Clarke
# 

import numpy as np
from latqcdtools.math.optimize import solve, persistentSolve
from latqcdtools.testing import print_results, concludeTest


def LHS1(x):
    return np.sin(x) - 1.


def testOptimize():

    lpass = True
    lpass *= print_results( solve(LHS1,guess=3,tol=1e-18,method='newton_krylov'), np.pi/2, prec=1e-7, text='newton_krylov' )
    lpass *= print_results( solve(LHS1,guess=3,tol=1e-8 ,method='fsolve')       , np.pi/2, prec=1e-7, text='fsolve'        )
    lpass *= print_results( persistentSolve(LHS1,guess=3,tol=1e-18)             , np.pi/2, prec=1e-7, text='persistent'    )

    concludeTest(lpass) 



if __name__ == '__main__':
    testOptimize()
