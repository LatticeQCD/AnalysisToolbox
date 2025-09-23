# 
# math.py                                                               
# 
# D. Clarke
# 
# Wrappers for some math functions, some math functions that I couldn't find in numpy or scipy, and some methods
# for comparing math objects. Some functions already existing in numpy or scipy may be implemented when I want
# to expose what those functions do under the hood.
#


import numpy as np
import scipy as sp
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import isArrayLike, cleanOutput
from latqcdtools.base.check import checkType


# ----------------------------------------------------------------- MATRICES


def id(N) -> np.ndarray:
    """ 
    NxN complex identity matrix
    """
    checkType('int',N=N)
    return np.eye(N,dtype=complex)


def ze(N) -> np.ndarray:
    """ 
    NxN complex zero matrix
    """
    checkType('int',N=N)
    return np.zeros((N,N), dtype=complex)


def isMatrix(mat) -> bool:
    checkType(np.ndarray,mat=mat)
    return mat.ndim==2


def checkMatrix(mat):
    if not isMatrix(mat):
        logger.TBRaise('Expected matrix. Got shape',np.shape(mat))


def isVector(vec) -> bool:
    checkType(np.ndarray,mat=vec)
    return vec.ndim==1


def checkVector(vec):
    if not isVector(vec):
        logger.TBRaise('Expected vector. Got shape',np.shape(vec))


def isSquare(mat) -> bool:
    checkMatrix(mat)
    return mat.shape[0] == mat.shape[1]


def checkSquare(mat):
    """ 
    Make sure mat is a square np.ndarray object. 

    Args:
        mat (np.ndarray)
    """
    if not isSquare(mat):
        logger.TBRaise('Expected square matrix. Got shape',np.shape(mat))


def dagger(arr) -> np.ndarray:
    checkType(np.ndarray,mat=arr)
    if arr.ndim > 2:
        logger.TBRaise('Expected ndim < 3. Got ndim =',arr.ndim) 
    return np.conjugate(arr).T


def isUnitary(mat) -> bool:
    """
    Unitary matrices U satisfy U^dag U = 1.
    
    Args:
        mat (np.ndarray)

    Returns:
        bool: True if unitary
    """
    checkMatrix(mat)
    N = len(mat[0])
    return rel_check(id(N),mat @ dagger(mat))


def isSpecial(mat) -> bool:
    """
    Special matrices M satisfy det(M) = 1.

    Args:
        mat (np.ndarray)

    Returns:
        bool: True if special
    """
    checkSquare(mat)
    return rel_check(np.linalg.det(mat), 1.)


def isSymmetric(mat) -> bool:
    """
    Symmetric matrices satisfy M^t = M.

    Args:
        mat (np.ndarray)

    Returns:
        bool: True if symmetric
    """
    checkSquare(mat)
    return np.allclose(mat, mat.T)


def isHermitian(mat) -> bool:
    """
    Hermitian matrices satisfy dagger(M)=M

    Args:
        mat (np.ndarray)

    Returns:
        bool: True if hermitian 
    """
    checkSquare(mat)
    return np.allclose(mat, dagger(mat))


def isAntihermitian(mat) -> bool:
    """
    Antihermitian matrices satisfy dagger(M)=-M

    Args:
        mat (np.ndarray)

    Returns:
        bool: True if antihermitian 
    """
    checkSquare(mat)
    return np.allclose(mat, -dagger(mat))


def isOrthogonal(mat) -> bool:
    """
    Orthogonal matrices satisfy M^t M = id

    Args:
        mat (np.ndarray)

    Returns:
        bool: True if orthogonal 
    """
    checkSquare(mat)
    N = len(mat[0])
    return np.allclose(mat@mat.T, id(N))


def isHankel(mat) -> bool:
    """
    Hankel matrices look like (3d example)
    a b c
    b c d
    c d e

    Args:
        mat np.ndarray 

    Returns:
        bool: True if Hankel 
    """
    checkSquare(mat)
    hankel = True
    for j in range(mat.shape[0]):
        for i in range(j):
            for k in range(j-i):
                hankel *= mat[i,j]==mat[i+k,j-k]
    return hankel


def TA(mat) -> np.ndarray:
    """
    Make mat traceless and antihermitian

    Args:
        mat (np.ndarray)

    Returns:
        bool: _description_
    """
    checkSquare(mat)
    return 0.5*(mat-dagger(mat))


def isPositiveSemidefinite(mat,details=False,eps=1e-12) -> bool:
    """ Returns true if mat is positive semidefinite. Otherwise, if details=True,
    list the eigenvalues that are not >=0.

    Args:
        mat (np.ndarray)
        details (bool, optional): If true, report problematic eigenvalues. Defaults to False.

    Returns:
        bool: True if positive semidefinite 
    """
    checkSquare(mat)
    eigenvalues = np.linalg.eigvals(mat)
    positiveSemidefinite = np.all(eigenvalues >= -eps)
    if details and (not positiveSemidefinite):
        logger.info('Problem eigenvalues:')
        for i in range(len(eigenvalues)):
            if not eigenvalues[i]>=0:
                logger.info(f'  lambda[{i}] = {eigenvalues[i]}')
    return positiveSemidefinite


def regulate(mat,svdcut=1e-12) -> np.ndarray:
    """ 
    If a matrix's singular values are too small, it will be ill-conditioned,
    making it difficult to invert and hence reducing numerical stability. This method
    extracts its singular values using SVD, then doctors the singular values to reduce
    the condition number. In the context of applying an SVD cut to a covariance 
    matrix, see e.g. Appendix D of 10.1103/PhysRevD.100.094508.

    Args:
        mat (np.ndarray)
        svdcut (float, optional): condition number threshold. Defaults to 1e-12.

    Returns:
        np.ndarray: regulated matrix 
    """
    logger.debug('regulate matrix with condition number',cleanOutput(np.linalg.cond(mat)))
    U, s, Vdagger = sp.linalg.svd(mat)
    smax = np.max(np.linalg.eigvals(mat))
    s[s < svdcut*smax] = svdcut*smax
    res = U @ np.diag(s) @ Vdagger 
    logger.debug('new condition number:',cleanOutput(np.linalg.cond(res)))
    return res


def invert(mat,method='scipy',svdcut=1e-12) -> np.ndarray:
    """ 
    Invert matrix.

    Args:
        mat (np.ndarray): to-be-inverted matrix 
        method (str): algorithm for inverting the matrix
        
    Returns:
        np.ndarray: mat^{-1} 
    """
    checkType(np.ndarray,mat=mat)
    checkSquare(mat)
    checkType(str,method=method)
    if method=='scipy':
        return sp.linalg.inv(mat)
    elif method=='numpy': # Seems to be quite slow
        return np.linalg.inv(mat)
    elif method=='pinv': 
        return np.linalg.pinv(mat)
    elif method=='svd':
        U, s, Vdagger = sp.linalg.svd(mat)
        return dagger(Vdagger) @ np.diag(1/s) @ dagger(U)
    elif method=='auto':
        condition_number = np.linalg.cond(mat)
        if condition_number > 1/svdcut:
            res = invert(regulate(mat,svdcut=svdcut),'pinv')
            return res
        else:
            return invert(mat,'scipy')
    else:
        logger.TBRaise('Unrecognized inverter',method)


def logDet(mat) -> float:
    """ 
    Logarithm of determinant. 
    """
    checkType(np.ndarray,mat=mat)
    checkSquare(mat)
    _, ans = np.linalg.slogdet(mat)
    return ans


def exp(mat) -> np.ndarray:
    """
    Exponential of square matrix.
    """
    checkSquare(mat) 
    return sp.linalg.expm(mat)


def log(mat) -> np.ndarray:
    """
    Natural logarithm of square matrix.
    """
    checkSquare(mat)
    return sp.linalg.logm(mat)


def pow(mat,power) -> np.ndarray:
    """
    Matrix power. 
    """
    checkSquare(mat)
    checkType('int',power=power)
    return np.linalg.matrix_power(mat, power)


class SUN(np.ndarray):

    """ 
    A member of the Lie group SU(N). Implemented as a subclass of the np.ndarray class. This gives us access already
    to all the nice features of np.ndarray and lets us leverage the speed of numpy.
        g.trace()
        g.det()
        g.dagger()
        g[i,j], which can be used to access and assign
        g + h
        g*h = g@h
        2*g
    """

    def __new__(cls, N=None, mat=None):
        if N is None:
            logger.TBRaise('Must specify N')
        if mat is None:
            mat = ze(N) 
        obj = np.asarray(mat, dtype=complex)
        if obj.shape != (N, N):
            logger.TBRaise(f"SU({N}) matrices must have shape ({N},{N})")
        return np.copy(obj).view(cls)


    def __repr__(self) -> str:
        return "SU(N)"


    def __mul__(self, other):
        # Perform matrix multiplication instead of element-wise multiplication
        return np.dot(self, other)


    def __pow__(self,power):
        # Perform matrix power instead of element-wise power
        return pow(self, power)


    def trace(self):
        return super().trace().item()


    def dagger(self):
        return dagger(self) 


    def det(self):
        return np.linalg.det(self)


    def isSUN(self) -> bool:
        """ 
        Check that I have det=1 and am unitary. 
        """
        if not (isSpecial(self) and isUnitary(self)):
            return False
        return True


    def setToMatrix(self,other):
        """ 
        Turn into RHS link. 
        """
        np.copyto(self,np.asarray(other,dtype=complex))


    def setToZero(self):
        """ 
        Turn into zero matrix. 
        """
        N = len(self)
        self.setToMatrix(ze(N))


    def setToIdentity(self):
        """ 
        Turn into identity matrix. 
        """
        N = len(self)
        self.setToMatrix(id(N)) 


# --------------------------------------------------------------- OTHER MATH


def pnorm(arr,p=2) -> float:
    """
    Returns p-norm of vector or matrix arr.

    Args:
        arr (np.ndarray)
        p (float, optional): Defaults to 2.

    Returns:
        float: pnorm 
    """
    checkType(np.ndarray,arr=arr)
    checkType("real",p=p)
    if p < 1:
        logger.TBRaise('p-norm defined only for p>=1.')
    if arr.ndim > 2: 
        logger.TBRaise('Expected ndim < 3. Got ndim =',arr.ndim) 
    if p==1:
        if isMatrix(arr):
            res = np.max(np.sum(np.abs(arr), axis=0))
        elif isVector(arr):
            res = np.sum(np.abs(arr))
    elif p==2:
        if isMatrix(arr):
            res = np.sqrt(np.trace(dagger(arr) @ arr)).real
        elif isVector(arr):
            res = np.sqrt(         dagger(arr) @ arr).real
    elif p==np.inf:
        if isMatrix(arr):
            res = np.max(np.sum(np.abs(arr), axis=1))
        elif isVector(arr):
            res = np.max(np.abs(arr))
    else:
        res = np.sum(np.abs(arr)**p)**(1./p) 
    return res


def normalize(arr,p=2) -> np.ndarray:
    """
    Normalize vector or matrix arr using p-norm.

    Args:
        arr (np.ndarray)
        p (float, optional): Defaults to 2.

    Returns:
        np.ndarray: normalized array 
    """
    checkType(np.ndarray,arr=arr)
    return arr/pnorm(arr,p)


def fallFactorial(n,m) -> float:
    """ 
    Falling factorial n fall to m. 
    """
    checkType('int',n=n)
    checkType('int',m=m)
    if m>n:
        logger.TBRaise("m>n.")
    return sp.special.poch(n-m+1,m)


def riseFactorial(n,m) -> float:
    """ 
    Rising factorial n rise to m. 
    """
    checkType('int',n=n)
    checkType('int',m=m)
    if n>m:
        logger.TBRaise("n>m.")
    return sp.special.poch(n,m)


def RMS(data) -> float:
    """
    Root-mean-square of data

    Args:
        data (np.ndarray)

    Returns:
        float: RMS of data 
    """
    checkType(np.ndarray,data=data)
    return np.sqrt(np.mean(data**2))


def quadrature(data) -> float:
    """
    Add data in quadrature

    Args:
        data (np.ndarray)

    Returns:
        float: data added in quadrature 
    """
    checkType(np.ndarray,data=data)
    return np.sqrt(np.sum(data**2))


def rel_check(a, b, prec = 1e-6, abs_prec = 1e-14) -> bool:
    """ 
    Check whether a and b are equal. a and b can be array-like, float-like, or complexes. If a
    and b are array-like, we check that they are element-wise equal within the tolerance. 

    Args:
        a (obj)
        b (obj)
        prec (float, optional): Relative precision. Defaults to 1e-6.
        abs_prec (float, optional): Absolute precision. Defaults to 1e-14.

    Returns:
        bool: True if a and b are equal. 
    """
    if isArrayLike(a):
        if np.shape(a) != np.shape(b):
            logger.TBRaise('a and b must have the same shape. Received a, b shapes =',np.shape(a),np.shape(b))
        return np.allclose( a, b, rtol = prec, atol = abs_prec)
    else:
        checkType("scalar",a=a)
        checkType("scalar",b=b)
        try:
            return np.isclose( a, b, rtol = prec, atol = abs_prec)
        except TypeError:
            logger.TBRaise('Expected reals, complexes, or array-like. Received a, b types =',type(a),',',type(b),exception=TypeError)
