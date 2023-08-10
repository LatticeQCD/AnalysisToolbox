# Bootstrapping routines

Given a set of $N$ measurements $\{x_1,...,x_N\}$, the statistical bootstrap allows one to estimate the error in
some function of the measurements $f$. Sometimes this is advantageous to error propagation, since analytically
calculating the error in the original function is too complicated. In the context of lattice
field theory, this happens e.g. when fitting correlators and trying to get the error from a fit parameter.
In the Analysistoolbox, these methods can be found in
```Python
import latqcdtools.statistics.bootstr
```
Just as with the [jackknife](jackknife.md), we stress that an advantage of the bootstrapping routines is that
you can pass them arbitrary functions.

## Ordinary bootstrap

Starting with our original measurements, one builds a bootstrap sample by drawing $N$ data from the original sample
with replacement. One repeats this process $K$ times. From bootstrap sample $i$, one gets an estimate of the mean
of interest. Averaging the $K$ means from each bootstrap sample gives a bootstrap mean.
The method
```Python
bootstr(func, data, numb_samples, sample_size = 0, same_rand_for_obs = False, conf_axis = 1, return_sample = False,
        seed = None, err_by_dist = False, args=(), nproc=DEFAULTTHREADS)
```
accomplishes this for an arbitrary $f$ `func`.

By default, the bootstrap sample size is equal to the original number of measurements. We resample with replacement.
The size can be adjusted with the `sample_size` argument.
By default the bootstrap is [parallelized](../base/speedify.md) with `DEFAULTTHREADS`
processes. Set `nproc=1` if you want to turn off parallelization.

## Gaussian bootstrap

The Gaussian bootstrap method,
```Python
bootstr_from_gauss(func, data, data_std_dev, numb_samples, sample_size = 1, same_rand_for_obs = False,
                   return_sample = False, seed = None, err_by_dist = True, useCovariance = False,
                   Covariance = None, args = (), nproc = DEFAULTTHREADS, asym_err=False)
```
will resample as follows: For each element of `data`, random data will be drawn from normal distributions
with means equal to the values in `data` and standard deviations from `data_std_dev`. This defines one
Gaussian bootstrap sample, and the function `func` is applied to the sample. This process is repeated
`numb_samples` times.

By default, the Gaussian bootstrap returns the median and 68-percentiles from the sample. You can return
the standard deviation instead by switching `err_by_dist` to `False`. You also have the option to get
back asymmetric quantiles/errors using `asymm_err=True`.
