latqcdtools.math.num_int
=============

```Python
integrateData(xdata, ydata, method='trapezoid'):
'''
Wrapper to integrate data. 

Args:
    xdata (array-like)
    ydata (array-like)
    method (str, optional): Integration method. Defaults to 'trapezoid'. Possibilities include:
                            > 'simpson' : Simpson rule.
                            > 'trapezoid' : Trapezoidal rule.

Returns:
    float: Area under ydata. 
'''
```
```Python
integrateFunction(func, a, b, method='persistent', args=(), stepsize=None, limit=1000, epsrel=1.49e-08, epsabs=1.49e-08, floatT=<class 'numpy.float64'>):
'''
Wrapper to integrate functions. Allows to conveniently adjust the stepsize, and can vectorize scipy.quad, 
which otherwise does not like to handle numpy arrays.

Args:
    func (func): Integrand. 
    a (float): Lower integration limit.
    b (float): Upper integration limit.
    method (str,optional): Integration method. Defaults to 'persistent_quad_trap'. Possibilities are:
                            > 'persistent' : Try various methods until something works. 
                            > 'quad' : Gaussian quadrature.
                            > 'trapezoid' : Trapezoidal rule.
    args (tuple, optional): Arguments to func. Defaults to ().
    stepsize (float, optional): _description_. Defaults to None.
    limit (int, optional): Upper bound on number of subintervals used in the adaptive algorithm. Defaults to 1000.
    epsrel (float, optional): Relative error tolerance. Defaults to 1.49e-8.
    epsabs (float, optional): Absolute error tolerance. Defaults to 1.49e-8.

Returns:
    float: Integral of func from a to b. 
'''
```
```Python
solveIVP(dydt, t0, tf, y0, method='RK45', args=(), epsrel=1.49e-08, epsabs=1.49e-08) -> numpy.ndarray:
'''
Wrapper to solve an initial value problem of the form

dy/dt = dydt(t, y)
y(t0) = y0

Args:
    dydt (func): RHS of IVP. Must have signature dydt(t,y).
    t0 (float): initial time 
    tf (float): final time 
    y0 (array-like): y(t0)
    method (str, optional): Integration method. Defaults to 'RK45'.
    args (tuple, optional): Arguments to dydt. Defaults to ().
    epsrel (float, optional): Relative error tolerance. Defaults to 1.49e-8.
    epsabs (float, optional): Absolute error tolerance. Defaults to 1.49e-8.

Returns:
    array-like: y(tf) 
'''
```
