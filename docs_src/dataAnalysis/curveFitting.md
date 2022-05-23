# Curve Fitting


There are many ways to fit a curve. There are strategies that minimize $\chi^2/{\rm d.o.f.}$ 
for when you know the functional form ahead of time, splines for when you don't, and other methods. 
The AnalysisToolbox includes some routines that are helpful for this purpose.

## $\chi^2$ minimization

In the module
```Python
import latqcdtools.statistics.fitting
``` 
one finds a `Fitter` class for carrying out fits. After constructing a `fitter` object, one can then use the associated 
methods to try various kinds of fits. These are generally wrappers from `scipy.optimize`. An easy example is given in 
`testing/fitting/simple_example.py`.

````{admonition} simple_example.py
 :class: toggle
```Python
#!/usr/bin/env python3
from latqcdtools.readin import *
from latqcdtools.fitting import *

print("\n Example of a simple 3-parameter quadratic fit.\n")

def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a * x**2 + b * x + c

xdata, ydata, edata = read_in("wurf.dat", 1, 3, 4)

fitter = Fitter(fit_func, xdata, ydata, expand = False)

res, res_err, chi_dof, pcov = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['levenberg', 'curve_fit'],
                                             ret_pcov = True)

print(" a , b,  c : ",res)
print(" ae, be, ce: ",res_err)
print("chi2/d.o.f.: ",chi_dof)
print("       pcov: \n",pcov,"\n")

fitter.plot_fit()
plt.show()
```
````

If one would rather use gnuplot to fit for some reason, there is also in `fitting.py` a python 
wrapper (`fit_gnuplot`) for passing fitting commands to gnuplot. An example of how to use this 
is given in `testing/fitting/test_fit_gnuplot`.

**IMPORTANT: The functions that you pass to these fitting routines have to be able to handle arrays!** 
E.g. you pass it `[x0, x1, ..., xN]` and get back `[f(x0), f(x1), ..., f(xN)]`. It is written this 
way to force better performance; if it were a typical loop it would be slow. If you are having 
trouble figuring out how to write your function in a way to handle arrays, a good starting point 
can be to use [np.vectorize](https://numpy.org/doc/stable/reference/generated/numpy.vectorize.html).

### Levenberg

The [Levenberg-Marquardt](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) 
algorithm seems to be pretty popular. It varies between Newton-Raphson and steepest descent. 
Normally one needs to know also the derivative of the to-be-fitted function. One can carry out 
such an "exact" Levenberg-Marquardt by passing the gradient and hessian fit functions to `do_fit`. 
However the AnalysisToolbox also gives you the option not to pass the gradient; in this case the 
gradient will be estimated numerically. Examples of the exact and numerical are respectively 
shown below.

````{admonition} Levenberg_example.py 
 :class: toggle
```Python
def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a * x**2 + b * x + c

def grad_fit_func(x, a, b, c):
    return [x**2, x, 1]

def hess_fit_func(x, a, b, c):
    return [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

xdata, ydata, edata = read_in("wurf.dat", 1, 3, 4)

# Exact Levenberg
fitter = Fitter(fit_func, xdata, ydata, edata, grad = grad_fit_func, hess = hess_fit_func,
                        func_sup_numpy=False)
res, res_err, chi_dof = fitter.do_fit(algorithm="levenberg", start_params = [1, 1, 1])

# Numerical Levenberg
fitter = fitting.Fitter(fit_func, xdata, ydata, edata, func_sup_numpy = True)
res, res_err, chi_dof = fitter.do_fit(algorithm="levenberg", start_params = [1, 1, 1])
```
````

### Other algorithms

`fitting.py` also has wrappers for several algorithms implemented in `scipy.optimize`. These include
- L-BFGS-B
- TNC
- Powell
- COBYLA
- SLSQP

If one does not specify an algorithm, the `try_fit` method will attempt all of them and return the 
result of whichever one had the best $\chi^2/{\rm d.o.f.}$

## Splines

There are several methods in the toolbox to fit a 1D spline to some `xdata` and `ydata`, possibly 
with y-error bars `edata`. These can be found in `latcqdtools/spline_interpolate.py`. One of the 
most straightforward methods is `bspline`, which can be called as, e.g.
```Python
xspline, yspline, yespline = bspline(xdata, ydata, edata, xmin=0.0 xmax=1.0, order= 4)
```
The outputs `xspline`, `yspline`, and `yespline` are `numpy` arrays that together make the graph 
of the spline along with its error bands. `xmin` and `xmax` specify the range over which the data 
should be fit, and `order` is the order of the spline.

Other options for fitting splines to data include `bspline_sci` and `LSQbspline_sci`, which are 
wrappers for the `scipy` method `UnivariateSpline` and `LSQUnivariateSpline`, respectively. However 
these two methods will not return an error band automatically, which means you need some way to 
assign errors yourself. If you need access to the spline itself, your best bet is to use 
`UnivariateSpline` and `LSQUnivariateSpline` directly, as these return the spline function.

### Fitting constrained splines

A simple example on how to fit a constrained spline (e.g. to some correlator data) is provided 
in the following:

````{admonition} spline_example.py 
 :class: toggle
```Python
#!/usr/bin/env python3
from latqcdtools import spline_interpolate
from latqcdtools import fitting
from latqcdtools import plotting
from matplotlib import pyplot

#Define the spline function with unknown coefficients. Constraint: at x=0.5 the 1st derivate shall be 0. The highest degree polynomial for order=4 is x^3.
myknots=[0.25] #Note: It is possible to let the fitter fit the knot positions by including them in the params array (see below, e.g. knots=params[0:2], coeffs=params[2:])
myorder=4
def mysplinefunc(x, params):
    return spline_interpolate.constraint_spline(x, knots=myknots, coeffs=params, base_point=0, order=myorder, constraints=[[0.5,1,0]])

#our correlator data
xdata = [0.0688, 0.1085, 0.1616, 0.2247, 0.2914, 0.3571, 0.4178, 0.4527]
ydata = [0.0098, 1.3194, 1.9954, 2.3360, 2.6893, 2.9958, 3.4421, 3.5319]

#The spline has (order+nknots) parameters. This simply sets each start parameter to 1.
mystart_params=[1 for i in range(len(myknots)+myorder)]

fitter = fitting.Fitter(func=mysplinefunc, xdata=xdata, ydata=ydata, expand=False)
fitparams, fitparams_err, chi_dof, cov = fitter.do_fit(start_params = mystart_params, ret_pcov = True) 

#Visualize the results
print(fitparams) 
print(fitparams_err)
print(chi_dof)
print(cov)

plotting.plot_dots(xdata, ydata)
plotting.plot_func(mysplinefunc, args = (fitparams,), xmax=0.5)

pyplot.show()
```
````
