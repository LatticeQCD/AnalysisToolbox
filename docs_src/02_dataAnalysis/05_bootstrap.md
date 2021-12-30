# Bootstrap 

Quantities such as masses are extracted by fitting some function plotted 
against an independent variable, here called $r$. Sometimes it is not clear what a good fitting range 
for $r$ is. For example Debye masses can be extracted by fitting to an exponential that is only valid 
for long distances, and there can be ambiguity in selecting an $r_{\text{min}}$. In such a case, one 
may obtain multiple estimates $m(r_{\text{min}})$ for different $r_{\text{min}}$ that are all similar 
to each other and highly correlated. The Gaussian bootstrap allows one to obtain an average and error 
bar under these circumstances. The method
```Python
avg, err = bootstr_add_dist(data, errors, nstat = 1000, plot_hist = False)
```
works as follows: Given possibly correlated `data` and corresponding error bars `errors` that are 
assumed to be Gaussian, resample by drawing for data point `i`, `nstat` new resampled measurements 
from a Gaussian distribution with mean `data[i]` and standard deviation `errors[i]`. Concatenate 
these `ndata*nstat` resampled measurements into a new distribution. The average `avg` is taken as 
the median of this new distribution, and the error `err` is the distance between the median and 
the 68% quantile.
