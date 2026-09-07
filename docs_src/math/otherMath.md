# Other mathematics 

## Solving equations numerically

In the module
```Python
latqcdtools.math.optimize
```
there is a the method
```Python
persistentSolve(LHS, guess, tol=1e-8, maxiter=200)
```
which will try to solve the equation `LHS==0` within tolerance `tol`, using up to 
`maxiter` iterations. This tries a few SciPy methods: in order,
`newton_krylov`, `fsolve`, then `root`. This is not necessarily the most optimal
order. It stops when one of them succeeds.

## Constructing polynomials

The module
```Python
latqcdtools.math.polynomials
```
contains `Polynomial` and `Rational` objects that can be used to succinctly represent
polynomials or rational functions. For example
```Python
p = Polynomial([A0, 0., A2, 0. A4])
p(x)
```
constructs a polynomial of only even powers up to fourth order.

## Padé approximants

The module
```Python
latqcdtools.math.pade
```
contains the `singlePointPade` class, which builds a rational approximation from the
Taylor coefficients of a function about a single point. Given coefficients `c` with
`f(x) = sum_k c[k] (x-x0)**k`, the call
```Python
R = singlePointPade(c, p=3, q=3, x0=0.)
R(x)
```
returns a callable rational function `P(x-x0)/Q(x-x0)` with `deg(P)=p`, `deg(Q)=q`,
whose own Taylor expansion about `x0` matches `f` through order `p+q`. The denominator is
normalized so `Q(0)=1`, and the object is backed by a `Rational` (see above). Both `p` and
`q` are required, and must satisfy `p+q <= len(c)-1`. Under the hood this wraps
`scipy.interpolate.pade`.

## Special functions

Most special functions are covered by SciPy, but some either somehow return extra values
or have notation that David is not used to. Therefore you can find
- `riseFactorial`: Compute $(n)^m$.
- `fallFactorial`: Compute $(n)_m$.
- `logDet`: Compute logarithm of determinant of a matrix.
