# Bootstrapping routines

Given a set of $N$ measurements $\{x_1,...,x_N\}$, the statistical bootstrap allows one to estimate the error in
some function of the measurements $f$. Sometimes this is advantageous to error propagation, since analytically
calculating the error in the original function is too complicated. In the context of lattice
field theory, this happens e.g. when fitting correlators and trying to get the error from a fit parameter.
In the Analysistoolbox, these methods can be found in
```Python
import latqcdtools.statistics.bootstr
```

## Ordinary bootstrap

Starting with our original measurements, one builds a bootstrap sample by drawing $N$ data from the original sample
with replacement. One repeats this process $K$ times. From bootstrap sample $i$, one gets an estimate of the mean
of interest. Averaging the $K$ means from each bootstrap sample gives a bootstrap mean.
The method
```Python
bootstr(func, data, numb_samples, sample_size = 0, same_rand_for_obs = False, conf_axis = 1, return_sample = False,
        seed = None, err_by_dist = False, args=(), nproc=DEFAULTTHREADS):
```
accomplishes this for an arbitrary $f$ `func`.

By default, the bootstrap sample size is equal to the original number of measurements. We resample with replacement.
The size can be adjusted with the `sample_size` argument.

## Gaussian bootstrap

Quantities such as masses are extracted by fitting some function plotted 
against an independent variable, here called $r$. Sometimes it is not clear what a good fitting range 
for $r$ is. For example Debye masses can be extracted by fitting to an exponential that is only valid 
for long distances, and there can be ambiguity in selecting an $r_{\text{min}}$. In such a case, one 
may obtain multiple estimates $m(r_{\text{min}})$ for different $r_{\text{min}}$ that are all similar 
to each other and highly correlated. The Gaussian bootstrap
allows one to obtain an average and error bar under these circumstances. The method
```Python
avg, err = bootstr_add_dist(data, errors, nstat = 1000, plot_hist = False)
```
works as follows: Given possibly correlated `data` and corresponding error bars `errors` that are 
assumed to be Gaussian, resample by drawing for data point `i`, `nstat` new resampled measurements 
from a Gaussian distribution with mean `data[i]` and standard deviation `errors[i]`. Concatenate 
these `ndata*nstat` resampled measurements into a new distribution. The average `avg` is taken as 
the median of this new distribution, and the error `err` is the distance between the median and 
the 68% quantile.
