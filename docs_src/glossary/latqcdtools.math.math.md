latqcdtools.math.math
=============

`RMS(data) -> float`

Root-mean-square of data

Args:
    data (np.ndarray)

Returns:
    float: RMS of data 

`TA(mat) -> numpy.ndarray`

Make mat traceless and antihermitian

Args:
    mat (np.ndarray)

Returns:
    bool: _description_

`checkMatrix(mat)`


`checkSquare(mat)`

Make sure mat is a square np.ndarray object. 

Args:
    mat (np.ndarray)

`dagger(mat) -> numpy.ndarray`


`fallFactorial(n, m) -> float`

Falling factorial n fall to m. 

`invert(mat, method='scipy', svdcut=1e-12) -> numpy.ndarray`

Invert matrix.

Args:
    mat (np.ndarray): to-be-inverted matrix 
    method (str): algorithm for inverting the matrix
    
Returns:
    np.ndarray: mat^{-1} 

`isAntihermitian(mat) -> bool`

Antihermitian matrices satisfy dagger(M)=-M

Args:
    mat (np.ndarray)

Returns:
    bool: True if antihermitian 

`isHankel(mat) -> bool`

Hankel matrices look like (3d example)
a b c
b c d
c d e

Args:
    mat np.ndarray 

Returns:
    bool: True if Hankel 

`isHermitian(mat) -> bool`

Hermitian matrices satisfy dagger(M)=M

Args:
    mat (np.ndarray)

Returns:
    bool: True if hermitian 

`isMatrix(mat) -> bool`


`isOrthogonal(mat) -> bool`

Orthogonal matrices satisfy M^t M = id

Args:
    mat (np.ndarray)

Returns:
    bool: True if orthogonal 

`isPositiveSemidefinite(mat, details=False, eps=1e-12) -> bool`
Returns true if mat is positive semidefinite. Otherwise, if details=True,
list the eigenvalues that are not >=0.

Args:
    mat (np.ndarray)
    details (bool, optional): If true, report problematic eigenvalues. Defaults to False.

Returns:
    bool: True if positive semidefinite 

`isSpecial(mat) -> bool`

Special matrices M satisfy det(M) = 1.

Args:
    mat (np.ndarray)

Returns:
    bool: True if special

`isSquare(mat) -> bool`


`isSymmetric(mat) -> bool`

Symmetric matrices satisfy M^t = M.

Args:
    mat (np.ndarray)

Returns:
    bool: True if symmetric

`isUnitary(mat) -> bool`

Unitary matrices U satisfy U^dag U = 1.

Args:
    mat (np.ndarray)

Returns:
    bool: True if unitary

`logDet(mat) -> float`

Logarithm of determinant. 

`normalize(arr)`


`quadrature(data) -> float`

Add data in quadrature

Args:
    data (np.ndarray)

Returns:
    float: data added in quadrature 

`regulate(mat, svdcut=1e-12) -> numpy.ndarray`

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

`rel_check(a, b, prec=1e-06, abs_prec=1e-14) -> bool`

Check whether a and b are equal. a and b can be array-like, float-like, or complexes. If a
and b are array-like, we check that they are element-wise equal within the tolerance. 

Args:
    a (obj)
    b (obj)
    prec (float, optional): Relative precision. Defaults to 1e-6.
    abs_prec (float, optional): Absolute precision. Defaults to 1e-14.

Returns:
    bool: True if a and b are equal. 

`riseFactorial(n, m) -> float`

Rising factorial n rise to m. 

