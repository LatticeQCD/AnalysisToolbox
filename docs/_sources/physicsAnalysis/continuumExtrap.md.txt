# Continuum-limit extrapolation 

Eventually we will have to take a continuum limit. The class for this is in
```Python
latqcdtools.physics.continuumExtrap
```
It is the `Extrapolator` class, which inherits from the `Fitter`,
described [here](../dataAnalysis/curveFitting.md). You instantiate an `Extrapolator` with
```Python
ext = Extrapolator(x, obs, obs_err, order=1, xtype="a", error_strat='propagation', nproc=DEFAULTTHREADS)
```
Here `x` either denotes lattice spacing data or $N_\tau$ data, which you can specify with the `xtype`.
We do power series extrapolations, with the `order` in $a^2$ specified by the user.
The `error_strat` is inherited from the `Fitter`; it decides whether you will compute errors
with error propagation or using the Hessian. This is [parallelized](../base/speedify.md) by
default, but you can set `nproc=1` to turn that off.

To perform an extrapolation, you can
```Python
ext.extrapolate(self,start_coeffs=None,prior=None,prior_err=None)
```
This gives you the option of passing priors, if you like. If you want to specify only some priors,
but not others, then pass `np.inf` as the corresponding `prior_err` for those you don't want to fix.
These will not be included in counting the degrees of freedom.

You can look at a plot with `ext.plot()`; you can pass it the same keyword arguments that
you can with the [plotter](../base/plotting.md). The class method `ext.showResults()` will
print the results to screen. 

For your convenience, this whole process is wrapped inside a method
```Python
continuumExtrapolate(x,obs,obs_err,order=1,show_results=False,plot_results=False,prior=None,
                     start_coeffs=None,prior_err=None,error_strat='propagation',xtype="a",
                     nproc=DEFAULTTHREADS):
```

