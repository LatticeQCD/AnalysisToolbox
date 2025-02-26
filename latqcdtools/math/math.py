# 
# math.py                                                               
# 
# D. Clarke
# 
# Wrappers for some math functions, some math functions that I couldn't find in numpy or scipy, and some methods
# for comparing math objects. 
#


import numpy as np
import scipy as sp
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import isArrayLike, cleanOutput 
from latqcdtools.base.check import checkType


id_2 = np.eye(2,dtype=complex)
id_3 = np.eye(3,dtype=complex)
id_4 = np.eye(4,dtype=complex)
ze_2 = np.zeros((3,3), dtype=complex)
ze_3 = np.zeros((3,3), dtype=complex)
ze_4 = np.zeros((3,3), dtype=complex)


def checkSquare(mat):
    """ 
    Make sure mat is a square np.ndarray object. 

    Args:
        mat (np.ndarray)
    """
    checkType(np.ndarray,mat=mat)
    if mat.shape[0] != mat.shape[1]:
        logger.TBRaise('Expected square matrix. Got shape',np.shape(mat))


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
        return np.conj(Vdagger.T) @ np.diag(1/s) @ np.conj(U.T)
    elif method=='auto':
        condition_number = np.linalg.cond(mat)
        if condition_number > 1/svdcut:
            res = invert(regulate(mat,svdcut=svdcut),'pinv')
            return res
        else:
            return invert(mat,'scipy')
    else:
        logger.TBRaise('Unrecognized inverter',method)


def isPositiveSemidefinite(mat,details=False,eps=1e-12) -> bool:
    """ Returns true if mat is positive semidefinite. Otherwise, if details=True,
    list the eigenvalues that are not >=0.

    Args:
        mat (np.ndarray)
        details (bool, optional): If true, report problematic eigenvalues. Defaults to False.

    Returns:
        bool: True if positive semidefinite 
    """
    checkType(np.ndarray,mat=mat)
    eigenvalues = np.linalg.eigvals(mat)
    positiveSemidefinite = np.all(eigenvalues >= -eps)
    if details and (not positiveSemidefinite):
        logger.info('Problem eigenvalues:')
        for i in range(len(eigenvalues)):
            if not eigenvalues[i]>=0:
                logger.info(f'  lambda[{i}] = {eigenvalues[i]}')
    return positiveSemidefinite 


def isSymmetric(mat) -> bool:
    checkType(np.ndarray,mat=mat)
    return np.allclose(mat, mat.T)


def normalize(arr):
    checkType(np.ndarray,arr=arr)
    return arr/np.sum(np.abs(arr))


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


def logDet(mat) -> float:
    """ 
    Logarithm of determinant. 
    """
    checkType(np.ndarray,mat=mat)
    checkSquare(mat)
    _, ans = np.linalg.slogdet(mat)
    return ans


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
            logger.TBRaise('Expected reals, complexes, or array-like. Received a, b types =',type(a),',',type(b))
