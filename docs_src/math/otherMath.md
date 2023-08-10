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

## Special functions

Most special functions are covered by SciPy, but some either somehow return extra values
or have notation that David is not used to. Therefore you can find
- `riseFactorial`: Compute $(n)^m$.
- `fallFactorial`: Compute $(n)_m$.
- `logDet`: Compute logarithm of determinant of a matrix.
