latqcdtools.statistics.bootstr
=============

```python
_autoSeed(seed) -> int:
"""
We use seed=None to flag the seed should be automatically chosen. The problem is that we need
seed to be an integer when enforcing that different bootstrap samples use different seeds. 
"""
```
```python
_biasCorrect(fbar, sampleval):
"""
Bias-corrected central value 2 f(xbar) - <f*>. Here f(xbar) is func of the full data, <f*> is the
average over bootstrap samples, and <f*> - f(xbar) estimates the O(1/N) bias of f(xbar). The result
carries a Monte Carlo error of about (bootstrap error)/sqrt(numb_samples), even for linear func. 
"""
```
```python
bootstr(func, data, numb_samples, same_rand_for_obs=False, conf_axis=1, return_sample=False, seed=None, err_by_dist=True, args=(), nproc=1):
"""
Bootstrap for arbitrary functions. This routine resamples the data and passes them to func in the same
format as the input, so func should compute an observable from a given data set. The central value is
the bias-corrected estimate 2 f(xbar) - <f*>, where f(xbar) is func of the full data and <f*> is the
average over bootstrap samples. Like the jackknife, this removes the O(1/N) bias of f(xbar). It carries
a Monte Carlo error of about (bootstrap error)/sqrt(numb_samples). The error is computed from the spread
of func over the bootstrap samples. The output of func may be a scalar or a numpy object. For multidimensional 
data, conf_axis says which axis holds the configurations to be resampled.

Args:
    func (callable): Function that calculates the observable.
    data (array-like): Input data.
    numb_samples (int): Number of bootstrap samples.
    same_rand_for_obs (bool, optional): Every index before conf_axis labels an observable. If True, all
      observables are resampled with the same random configurations, i.e. they are treated as perfectly
      correlated. If False, each observable gets its own random configurations. Indices after conf_axis
      always share the random configuration. Defaults to False.
    conf_axis (int, optional): Axis to resample. Defaults to 1; it is set to 0 for one-dimensional data.
    return_sample (bool, optional): Also return the results from the individual samples? Defaults to False.
    seed (int, optional): Seed for the random generator. Defaults to None, which chooses a random seed.
    err_by_dist (bool, optional): Take the error from the 68% quantiles instead of the standard deviation?
      The error is then the distance from the median of the bootstrap distribution to its quantiles, i.e.
      a width of the distribution, applied around the bias-corrected central value. Defaults to True.
    args (tuple or dict, optional): Extra arguments for func. A dict is passed as **args. Defaults to ().
    nproc (int, optional): Number of threads. Defaults to 1, i.e. no parallelization. Starting the
      parallel pool has an overhead, so increasing nproc only helps when func is slow, for example
      if you do a curve fit in each sample. Feel free to benchmark it for your use-case.

Returns:
    samples (optionally), bias-corrected central value, bootstrap error
"""
```
```python
bootstr_from_gauss(func, data, data_std_dev, numb_samples, same_rand_for_obs=False, return_sample=False, seed=None, err_by_dist=True, useCovariance=False, Covariance=None, args=(), nproc=1, asym_err=False):
"""
Gaussian (parametric) bootstrap. Like bootstr, but each sample is drawn from a normal distribution around
the mean values in data, with width data_std_dev or covariance Covariance. The central value is the
bias-corrected estimate 2 f(xbar) - <f*>, as in bootstr, with f(xbar) being func of the mean values.

Args:
    func (callable): Function that calculates the observable.
    data (array-like): Mean value of each observable.
    data_std_dev (array-like): Error of each value in data, e.g. the error of a mean.
    numb_samples (int): Number of bootstrap samples.
    same_rand_for_obs (bool, optional): Use the same random numbers for every observable, i.e. treat them
      as perfectly correlated? Defaults to False.
    return_sample (bool, optional): Also return the results from the individual samples? Defaults to False.
    seed (int, optional): Seed for the random generator. Defaults to None, which chooses a random seed.
    err_by_dist (bool, optional): Take the error from the 68% quantiles instead of the standard deviation?
      The error is then the distance from the median of the bootstrap distribution to its quantiles, i.e.
      a width of the distribution, applied around the bias-corrected central value. Defaults to True.
    useCovariance (bool, optional): Draw each observable from a multivariate normal? Defaults to False.
    Covariance (array-like, optional): Covariance matrix of each observable, used if useCovariance. Defaults
      to None, which uses diag(data_std_dev**2).
    args (tuple or dict, optional): Extra arguments for func. A dict is passed as **args. Defaults to ().
    nproc (int, optional): Number of threads. Defaults to 1, i.e. no parallelization. Starting the
      parallel pool has an overhead, so increasing nproc only helps when func is slow, for example
      if you do a curve fit in each sample. Feel free to benchmark it for your use-case.
    asym_err (bool, optional): Return both distances from the median to the 68% quantiles, if err_by_dist,
      instead of the larger one. Defaults to False.

Returns:
    samples (optionally), bias-corrected central value, bootstrap error
"""
```
```python
class nimbleBoot(func, data, numb_samples, same_rand_for_obs, conf_axis, return_sample, seed, err_by_dist, args, nproc):
```
```python
class nimbleGaussianBoot(func, data, data_std_dev, numb_samples, same_rand_for_obs, return_sample, seed, err_by_dist, useCovariance, Covariance, args, nproc, asym_err):
```
