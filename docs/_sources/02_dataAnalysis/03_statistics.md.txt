# Statistics 

## Mean, median, and error 
There are wrappers for `np.mean` and `np.median` called `std_mean` 
and `std_median`. The advantage of using the wrappers is that you don't have to specify the Axis. 
There is also a `std_err`, which calculates the error bar of the mean of the array, assuming 
the data are independent.

## Gaussian difference test 
A Gaussian difference test can be used to check whether two 
measurements are statistically consistent with one another. Given are two measurements `x1` 
and `x2` drawn from Gaussian distributions with the same mean and respective error bars `e1` 
and `e2`. A call to
```Python
gaudif(x1,e1,x2,e2)
```
returns the q-value, which is the likelihood that `x1` and `x2` are at least as far apart as 
was observed. 

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
