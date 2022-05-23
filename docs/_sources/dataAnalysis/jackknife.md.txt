# Jackknife 

To generate jackknife error bars, please use the module
```Python
import latqcdtools.statistics.jackknife
```
The central method of this file is the `jackknife` method. A call to
```Python
jackknife(func, data, numb_blocks=20, conf_axis=1, args=(), cov=False)
```
returns a jackknife average, error, and (if desired) the covariance, in that order. Here `data` is 
a time series of raw data you are interested in, and `func` is some function of that data. 
`numb_blocks` is the number of jackknife blocks, `args` are the arguments to pass to `func`. Set 
`cov=True` if you want to get the covariance. The `conf_axis` is needed when one wants to pass a 
2D array as `data`, which can happen, for instance, if `func` depends on multiple observables. 
An example will be given below.

This jackknife is quite general. Something easy one might do with the jackknife is calculate the 
Polyakov loop susceptibility from Polyakov loops. A rule of thumb to use when using the jackknife 
code is to "write your function the way it appears using math symbols". For example with the 
susceptibility one calculates $\chi_P=N_\sigma^3\left(\langle P^2\rangle - \langle P\rangle\right)$, 
so one correspondingly defines
```Python
Ns=32
def susc(P):
  return Ns**3*( np.mean(P**2) - np.mean(P)**2 )
```
where $P$ is the time series of measurements of $P$, and then calls the jackknife by
```Python
Pmean, Perr = jackknife(susc, P)
```
Sometimes your function may depend on multiple observables; for instance one may be interested in 
the susceptibility of $|P|$, for which one needs the real and imaginary parts of $P$. One can 
accomplish this using 32 jackknife blocks with, for instance,
```Python
def suscA(data):
  ReP = data[0]
  ImP = data[1]
  return Ns**3*( np.mean(ReP**2+ImP**2) - np.mean(ReP)**2 )

AbsPmean, AbsPerr = jackknife(suscA, [ReP, ImP], 32, 1)
```
Note that these Polyakov loop examples are meant to be instructional. A complete set of functions 
measuring Polyakov loop observables is given in `polyakovTools.py`, described in part in the part
of the documentation for [physics](../physicsAnalysis/physicsAnalysis.md) modules. 
**WARNING:** Although the `jackknife` 
method is very general, one thing that cannot be done is passing a lambda function. This is because 
the `jackknife` is parallelized using `concurrent.futures`, which is not able to pickle 
lambda functions.
