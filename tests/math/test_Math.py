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
    isSpecial, isSymmetric, isHermitian, isAntihermitian, isOrthogonal, isHankel, TA, pnorm, \
    normalize, exp, pow, log, ze
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
              [ -0.28004516+0.098896j  ,  -0.35911944-0.71267198j,  1j                    ] ])

O = np.array([[ np.cos(1.23), np.sin(1.23)],
              [-np.sin(1.23), np.cos(1.23)]])

hank = np.array([[1,2,3,4,5],
                 [2,3,4,5,6],
                 [3,4,5,6,7],
                 [4,5,6,7,8],
                 [5,6,7,8,9]])

genmat = np.array([[1j,1j         ,np.pi+1j],
                   [1 ,2+np.pi*1j ,3       ]])

genvec = np.array([1,1+1j,np.pi])


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

    lpass *= print_results( pnorm(genmat,2)     , np.linalg.norm(genmat,'fro') , text='2-norm mat' )
    lpass *= print_results( pnorm(genvec,2)     , np.linalg.norm(genvec,None ) , text='2-norm vec' )
    lpass *= print_results( pnorm(genmat,1)     , np.linalg.norm(genmat,1)     , text='1-norm mat' )
    lpass *= print_results( pnorm(genvec,1)     , np.linalg.norm(genvec,1)     , text='1-norm vec' )
    lpass *= print_results( pnorm(genvec,np.inf), np.linalg.norm(genvec,np.inf), text='inf-norm vec' )
    lpass *= print_results( pnorm(genmat,np.inf), np.linalg.norm(genmat,np.inf), text='inf-norm mat' )

    test = ze(3) 
    for i in range(30):
        test += (1/math.factorial(i))*pow(A,i)
    lpass *= print_results( test, exp(A)     , text='matrix exponential' )
    lpass *= print_results( A   , log(exp(A)), text='matrix logarithm' )

    for p in [1,2,3,np.inf]:
        lpass *= print_results( pnorm(normalize(genvec,p),p), 1, text=f'{p} vec normalize')

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