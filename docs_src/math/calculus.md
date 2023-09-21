# Numerical calculus 

In the AnalysisToolbox we have some methods for numerical differentiation and integration.

## Numerical differentiation

Differentiation is implemented in
```Python
latqcdtools.math.num_deriv
```
We use central differences to calculate the numerical derivative. The method `best_h`
finds a finite difference that is small enough to suppress higher order differences,
but not so small that you suffer from errors due to machine precision.

To calculate $f'$ of $f$ you can call
```Python
diff_deriv(x, func, args = (), h = None)
```
We also have methods to compute the gradient, `diff_grad`, and the Hessian, `diff_hess`.

## Numerical integration

SciPy already has some pretty fast integrators, so we just wrap those. The wrappers can 
be found in
```Python
latqcdtools.math.num_int
```

One of the central wrappers is
```Python
integrateFunction(func,a,b,method='persistent',args=(),stepsize=None,limit=1000,epsrel=1.49e-8,epsabs=1.49e-8)
```
which carries out the one-dimensional integral of `func` from `a` to `b`. There are a variety of integration
methods to choose from that can be specified with `method`, including the trapezoid rule, the Romberg method,
and Gaussian quadrature. The solving continues until subsequent guesses fall within relative or absolute
differences `epsrel` and `epsabs`, or if it reaches `limit` evaluations. The default method `persistent`
tries various methods until something works.

Another wrapper is
```Python
integrateData(xdata,ydata,method='trapezoid')
```
which will find the area under `ydata` given the grid points `xdata`. This uses the trapezoid rule by default,
but it can also use the Simpson rule or the homebrew `spline` method, which fits the data with a spline,
then evaluates the area by quadrature.
