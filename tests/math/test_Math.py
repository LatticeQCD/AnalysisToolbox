# 
# testMath.py                                                               
# 
# D. Clarke
# 
# For testing some of the basic math functions.
#
import numpy as np
import math
import latqcdtools.base.logger as logger
from latqcdtools.math.math import fallFactorial, invert, RMS, isMatrix, isSquare, isUnitary, \
    isSpecial, isSymmetric, isHermitian, isAntihermitian, isOrthogonal, isHankel, TA
from latqcdtools.testing import print_results, concludeTest


logger.set_log_level('INFO')

mat = np.array( 
      [[ 1,2,-1],
       [ 2,1, 2],
       [-1,2, 1]])

tensor = np.ones((2,2,2))

U = np.array([[1 ,1  ],
              [1j,-1j]])/np.sqrt(2)

H = np.array([[ 0.17036869            ,  0.03794601-0.13535497j,  0.28004516-0.098896j  ],
              [ 0.03794601+0.13535497j,  0.01763769            ,  0.35911944+0.71267198j],
              [ 0.28004516+0.098896j  ,  0.35911944-0.71267198j, -0.8191237             ] ])

A = np.array([[  1j                    ,   0.03794601+0.13535497j,  0.28004516+0.098896j  ],
              [ -0.03794601+0.13535497j,   1j                    ,  0.35911944-0.71267198j],
              [ -0.28004516+0.098896j  ,  -0.35911944-0.71267198j,  1j             ] ])

O = np.array([[ np.cos(1.23), np.sin(1.23)],
              [-np.sin(1.23), np.cos(1.23)]])

hank = np.array([[1,2,3,4,5],
                 [2,3,4,5,6],
                 [3,4,5,6,7],
                 [4,5,6,7,8],
                 [5,6,7,8,9]])


def testMath():

    lpass = True

    lpass *= isMatrix(mat)
    lpass *= not isMatrix(tensor)
    lpass *= isSquare(mat)
    lpass *= not isSpecial(mat)
    lpass *= isSymmetric(mat)
    lpass *= isUnitary(U)
    lpass *= isHermitian(H)
    lpass *= isAntihermitian(A)
    lpass *= isOrthogonal(O)
    lpass *= isHankel(hank)

    testTA = TA(H)
    lpass *= isAntihermitian(testTA)
    lpass *= np.isclose(np.trace(testTA),0)

    lpass *= print_results(fallFactorial(23,23),1.*math.factorial(23),text='fallFactorial 23!')
    lpass *= print_results(fallFactorial(6,3),6*5*4,text='fallFactorial 6_fall_3')

    inv = invert(mat,'scipy')

    lpass *= print_results(inv,invert(mat,'numpy'),text='scipy vs numpy')
    lpass *= print_results(inv,invert(mat,'svd'),text='scipy vs svd')
    lpass *= print_results(inv,invert(mat,'pinv'),text='scipy vs pinv')

    data = np.array([7,12,48,1/2,np.pi])
    lpass *= print_results(22.392496977340823,RMS(data),text='RMS')

    concludeTest(lpass)


if __name__ == '__main__':
    testMath()