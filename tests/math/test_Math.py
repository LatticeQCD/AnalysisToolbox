# 
# testMath.py                                                               
# 
# D. Clarke
# 
# For testing some of the basic math functions.
#
import numpy as np
import math
from latqcdtools.math.math import fallFactorial, invert, RMS, isMatrix, isSquare, isUnitary, \
    isSpecial, isSymmetric, isHermitian, isAntihermitian, isOrthogonal, isHankel, TA, pnorm, \
    normalize, exp, pow, log, ze, regulate, eigh, applySVDCut
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger


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

    # Indefinite matrix: the cut is relative to the largest singular value 10, not to the largest eigenvalue 1e-3,
    # so the condition number is capped at 1/svdcut.
    lpass *= print_results(1e12,np.linalg.cond(regulate(np.diag([-10.,1e-3,1e-15]),svdcut=1e-12)),
                           text='regulate indefinite',prec=1e-10)

    # Symmetry checks are relative to the size of the matrix, so they do not depend on its units.
    asym = np.array([[1.,2.],[2.1,1.]])
    if isSymmetric(1e-20*asym) or not isSymmetric(1e-20*(asym+asym.T)):
        logger.TBFail('isSymmetric depends on scale')
        lpass = False

    # eigh returns orthonormal eigenvectors also for degenerate eigenvalues.
    J = np.ones((3,3))
    vals, vecs = eigh(1.6*np.eye(3)-0.6*J)
    lpass *= print_results(np.array([-0.2,1.6,1.6]),vals,text='eigh eigenvalues',prec=1e-14)
    lpass *= print_results(np.eye(3),vecs.T @ vecs,text='eigh orthonormal',prec=1e-14)

    lpass *= print_results(np.array([1e-6,1e-6,1.]),applySVDCut(np.array([1e-9,1e-6,1.]),1e-6),text='applySVDCut',prec=1e-14)

    # For an ill-conditioned matrix, invert 'auto' is the inverse of the regulated matrix.
    ill = np.array([[1.,1.],[1.,1.+1e-10]])
    lpass *= print_results(invert(regulate(ill,svdcut=1e-6),'scipy'),invert(ill,'auto',svdcut=1e-6),
                           text='invert auto',prec=1e-8)

    concludeTest(lpass)


if __name__ == '__main__':
    testMath()
