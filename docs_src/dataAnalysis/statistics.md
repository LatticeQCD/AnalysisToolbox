# Statistics 

A collection of useful methods for statistics calculations can be found in
```Python
import latqcdtools.statistics.statistics
```

## Mean, median, and standard error 
There are wrappers for `np.mean` and `np.median` called `std_mean` 
and `std_median`. The advantage of using the wrappers is that you don't have to specify the Axis. 
There is also a `std_err`, which calculates the error bar of the mean of the array, assuming 
the data are independent.

## Error propagation
The function  
```Python
error_prop_func(x, func, means, errors, grad=None, args=())
```
can be used to automatically propagate the `errors` of measurements `means` into the function `func`. Here `x`
is a variable the function depends on that does not have error. If you like, you can specify the gradient `grad`
of `func` yourself; otherwise this will be calculated numerically.

## Gaussian difference test (Z-test) 
A Gaussian difference test can be used to check whether two 
measurements are statistically consistent with one another. Given are two measurements `x1` 
and `x2` drawn from Gaussian distributions with the same mean and respective error bars `e1` 
and `e2`. A call to
```Python
gaudif(x1,e1,x2,e2)
```
returns the q-value, which is the likelihood that `x1` and `x2` are at least as far apart as 
was observed.

## Student different test (T-test)
In the case you have a small number of data, a more correct measure of the tension between two means
is given by the Student difference test. This can be called with
```Python
studif(x1,e1,ndat1,x2,e2,ndat2)
```
where `ndat1` is the number of measurements leading to mean `x1` and `ndat2` is the number of
measurements leading to mean `x2`. For large enough `ndat1` and `ndat2`, the results should
be similar to the Gaussian difference test.

## Integrated autocorrelation time
The integrated autocorrelation time is the ratio between 
the estimated variance of the sample mean and what this variance would have been if the data were 
independent; i.e. $\sigma^2_{\bar{X}}=\frac{\sigma^2}{N}\tau_{\text{int}}.$ 
A call to
```Python
getTauInt(timeSeries, nbins, tpickMax, acoutfileName)
```
returns an estimate for $\tau_{\text{int}}$, its error bar, its bias, and the Monte Carlo time 
separation at which it found this $\tau_{\text{int}}$. It takes a time series `timeSeries`, a 
number of jackknife bins needed for the calculation `nbins`, and an estimate for the largest that 
the autocorrelation time could be. (You can play around with this latter variable a bit if you 
are having trouble to get it to run. Usually I pick a little less than half the size of the series.) 
The results are saved by default in `acoutfileName=acor.d`, but you can change this as well.
