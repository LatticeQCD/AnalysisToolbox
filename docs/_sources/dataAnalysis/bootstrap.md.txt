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
of interest. The spread of these $K$ estimates gives the error. For the central value, note that $f$ of the
full sample, $f(\bar{x})$, has a bias of $\mathcal{O}(1/N)$ when $f$ is nonlinear, and the average over bootstrap
samples, $\langle f^*\rangle$, has about twice that bias. Their difference therefore estimates the bias, and the
routines return the bias-corrected value $2f(\bar{x})-\langle f^*\rangle$, just as the [jackknife](jackknife.md)
returns its bias-corrected mean. This value carries a Monte Carlo error of about $\sigma/\sqrt{K}$, where
$\sigma$ is the bootstrap error, even when $f$ is linear.
The method
```Python
bootstr(func, data, numb_samples, same_rand_for_obs = False, conf_axis = 1, return_sample = False,
        seed = None, err_by_dist = True, args=(), nproc=1)
```
accomplishes this for an arbitrary $f$ `func`.

Each bootstrap sample has the same size as the original number of measurements, and we resample with replacement.
The bootstrap can be [parallelized](../base/speedify.md) with the `nproc` argument. It is not parallelized by
default, since starting the parallel pool has an overhead that only pays off when `func` is slow, e.g. if
`func` carries out a fit.

## Gaussian bootstrap

The Gaussian bootstrap method,
```Python
bootstr_from_gauss(func, data, data_std_dev, numb_samples, same_rand_for_obs = False,
                   return_sample = False, seed = None, err_by_dist = True, useCovariance = False,
                   Covariance = None, args = (), nproc = 1, asym_err=False)
```
will resample as follows: For each element of `data`, random data will be drawn from normal distributions
with means equal to the values in `data` and standard deviations from `data_std_dev`. This defines one
Gaussian bootstrap sample, and the function `func` is applied to the sample. This process is repeated
`numb_samples` times.

The central value is bias-corrected as above, with $f(\bar{x})$ being $f$ of the values in `data`. By default, both bootstraps get the error from the
68-percentiles of the bootstrap distribution. You can use the standard deviation instead by switching
`err_by_dist` to `False`. You also have the option to get
back asymmetric quantiles/errors using `asym_err=True` (Gaussian bootstrap only).
