latqcdtools.math.num_deriv
=============

```python
_best_h(x):
"""
This routine attempts to automatically determine the best possible step size h when calculating numerical
derivatives. One does not like the step size to be too large, since one pays a ~h^2 penalty. On the other hand,
if the step size is too small, the function may not change within machine precision eps, and hence the derivative
will be incorrectly estimated.

Let M_0, M_3 > 0 s.t. |f| < M_0 and |f^(3)| < M_3. Let f_computed = f + f*eps, where eps is the machine precision.
For central differences,
    ( f(x+h) - f(x-h) )/2h - f' = h^3 f^(3)(x)/3.
Using the triangle inequality, one can construct an upper bound of the LHS (when calculated on a computer). The h
which minimizes this upper bound is
    h = (3 eps M_0/M_3)^(1/3).
Sadly this requires knowledge of M_3, which we presumably don't have since we resorted to numerical derivatives in
the first place. A compromise that has some memory of this optimization problem is
    h ~ x eps^(1/3).

We choose eps = 1.1e-16, which is roughly the difference between 1.0 and the next-smallest representable float less
than 1.0 in 64-bit. The difference between 1.0 and the next-largest float is slightly bigger. 
"""
```
```python
diff_deriv(x, func, args=(), h=None):
"""
Numerical derivative using central difference. 
"""
```
```python
diff_fit_grad(x, params, func, args=(), h=None):
"""
When fitting we're trying to optimize params, and hence we want to think of func as a function of
its parameters rather than x. 
"""
```
```python
diff_fit_hess(x, params, func, args=(), h=None):
"""
When fitting we're trying to optimize params, and hence we want to think of func as a function of
its parameters rather than x. 
"""
```
```python
diff_grad(params, func, args=(), h=None, floatT=numpy.float64) -> numpy.ndarray:
"""
Gradient using difference quotient. h may be None, a scalar, or an array with one step size per parameter.
"""
```
```python
diff_hess(params, func, args=(), h=None, floatT=numpy.float64) -> numpy.ndarray:
"""
Hessian using difference quotient. 
"""
```
```python
diff_jac(params, func, args=(), h=None) -> numpy.ndarray
```
```python
relative_h(x, scale) -> numpy.ndarray:
"""
Step sizes for each component of x that are relative to that component, so that a tiny x_i is never stepped
across zero. This matters for functions that are singular at zero, like 1/x. A zero x_i is stepped relative
to scale_i instead, e.g. its error bar. If that is also zero, fall back to _best_h(0). 

Args:
    x (np.ndarray): point at which to differentiate
    scale (np.ndarray): natural scale of each component, used when x_i = 0

Returns:
    np.ndarray: step sizes 
"""
```
