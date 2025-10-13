latqcdtools.math.math
=============

```Python
RMS(data) -> float:
'''
Root-mean-square of data

Args:
    data (np.ndarray)

Returns:
    float: RMS of data 
'''
```
```Python
TA(mat) -> numpy.ndarray:
'''
Make mat traceless and antihermitian

Args:
    mat (np.ndarray)

Returns:
    bool: _description_
'''
```
```Python
checkMatrix(mat)
```
```Python
checkSquare(mat):
'''
Make sure mat is a square np.ndarray object. 

Args:
    mat (np.ndarray)
'''
```
```Python
checkVector(vec)
```
```Python
dagger(arr) -> numpy.ndarray
```
```Python
exp(mat) -> numpy.ndarray:
'''
Exponential of square matrix.
'''
```
```Python
fallFactorial(n, m) -> float:
'''
Falling factorial n fall to m. 
'''
```
```Python
id(N) -> numpy.ndarray:
'''
NxN complex identity matrix
'''
```
```Python
invert(mat, method='scipy', svdcut=1e-12) -> numpy.ndarray:
'''
Invert matrix.

Args:
    mat (np.ndarray): to-be-inverted matrix 
    method (str): algorithm for inverting the matrix
    
Returns:
    np.ndarray: mat^{-1} 
'''
```
```Python
isAntihermitian(mat) -> bool:
'''
Antihermitian matrices satisfy dagger(M)=-M

Args:
    mat (np.ndarray)

Returns:
    bool: True if antihermitian 
'''
```
```Python
isHankel(mat) -> bool:
'''
Hankel matrices look like (3d example)
a b c
b c d
c d e

Args:
    mat np.ndarray 

Returns:
    bool: True if Hankel 
'''
```
```Python
isHermitian(mat) -> bool:
'''
Hermitian matrices satisfy dagger(M)=M

Args:
    mat (np.ndarray)

Returns:
    bool: True if hermitian 
'''
```
```Python
isMatrix(mat) -> bool
```
```Python
isOrthogonal(mat) -> bool:
'''
Orthogonal matrices satisfy M^t M = id

Args:
    mat (np.ndarray)

Returns:
    bool: True if orthogonal 
'''
```
```Python
isPositiveSemidefinite(mat, details=False, eps=1e-12) -> bool:
'''Returns true if mat is positive semidefinite. Otherwise, if details=True,
list the eigenvalues that are not >=0.

Args:
    mat (np.ndarray)
    details (bool, optional): If true, report problematic eigenvalues. Defaults to False.

Returns:
    bool: True if positive semidefinite 
'''
```
```Python
isSpecial(mat) -> bool:
'''
Special matrices M satisfy det(M) = 1.

Args:
    mat (np.ndarray)

Returns:
    bool: True if special
'''
```
```Python
isSquare(mat) -> bool
```
```Python
isSymmetric(mat) -> bool:
'''
Symmetric matrices satisfy M^t = M.

Args:
    mat (np.ndarray)

Returns:
    bool: True if symmetric
'''
```
```Python
isUnitary(mat) -> bool:
'''
Unitary matrices U satisfy U^dag U = 1.

Args:
    mat (np.ndarray)

Returns:
    bool: True if unitary
'''
```
```Python
isVector(vec) -> bool
```
```Python
log(mat) -> numpy.ndarray:
'''
Natural logarithm of square matrix.
'''
```
```Python
logDet(mat) -> float:
'''
Logarithm of determinant. 
'''
```
```Python
normalize(arr, p=2) -> numpy.ndarray:
'''
Normalize vector or matrix arr using p-norm.

Args:
    arr (np.ndarray)
    p (float, optional): Defaults to 2.

Returns:
    np.ndarray: normalized array 
'''
```
```Python
pnorm(arr, p=2) -> float:
'''
Returns p-norm of vector or matrix arr.

Args:
    arr (np.ndarray)
    p (float, optional): Defaults to 2.

Returns:
    float: pnorm 
'''
```
```Python
pow(mat, power) -> numpy.ndarray:
'''
Matrix power. 
'''
```
```Python
quadrature(data) -> float:
'''
Add data in quadrature

Args:
    data (np.ndarray)

Returns:
    float: data added in quadrature 
'''
```
```Python
regulate(mat, svdcut=1e-12) -> numpy.ndarray:
'''
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
'''
```
```Python
rel_check(a, b, prec=1e-06, abs_prec=1e-14) -> bool:
'''
Check whether a and b are equal. a and b can be array-like, float-like, or complexes. If a
and b are array-like, we check that they are element-wise equal within the tolerance. 

Args:
    a (obj)
    b (obj)
    prec (float, optional): Relative precision. Defaults to 1e-6.
    abs_prec (float, optional): Absolute precision. Defaults to 1e-14.

Returns:
    bool: True if a and b are equal. 
'''
```
```Python
riseFactorial(n, m) -> float:
'''
Rising factorial n rise to m. 
'''
```
```Python
ze(N) -> numpy.ndarray:
'''
NxN complex zero matrix
'''
```
```Python
class SUN(N=None, mat=None):
'''
A member of the Lie group SU(N). Implemented as a subclass of the np.ndarray class. This gives us access already
to all the nice features of np.ndarray and lets us leverage the speed of numpy.
    g.trace()
    g.det()
    g.dagger()
    g[i,j], which can be used to access and assign
    g + h
    g*h = g@h
    2*g
'''
```
