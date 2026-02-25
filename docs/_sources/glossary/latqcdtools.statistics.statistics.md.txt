latqcdtools.statistics.statistics
=============

```Python
AIC(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
'''
The Akaike information criterion (AIC) is a measure of how well a fit performs. It builds on the likelihood
function by including a penalty for each d.o.f. This is useful in a context where you have multiple models to
choose from,and hence different numbers of d.o.f. possible. It's also useful when you are worried about
overfitting. The preferred model minimizes the AIC.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): model parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
    
    Returns:
        float: AIC
'''
```
```Python
AICc(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
'''
Corrected AIC (AICc). When the sample size is smaller, it increases the chance AIC will select a model with too
many parameters. The AICc tries to further correct for this. In the limit that the number of data points goes to
infinity, one recovers the AIC.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): model parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
    
    Returns:
        float: corrected AIC
'''
```
```Python
BAIC(xdata, ydata, cov, func, args=(), params=(), Ncut=0, modelPrior=1) -> float:
'''
Bayesian Akaike information criterion of 2208.14983. It uses the chi^2 as its likelihood
function and includes penalties for having many fit parameters and cutting many data from
your original sample.

Args:
    xdata (array-like)
    ydata (array-like)
    cov (array-like): covariance matrix 
    func (func)
    args (tuple, optional): arguments to func. Defaults to ().
    params (tuple, optional): model parameters. Defaults to ().
    prior (array-like, optional): Bayesian priors. Defaults to None.
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
    Ncut (int, optional): The number of data trimmed from your fit. Defaults to 0.
    modelPrior (float, optional): The prior probability of your model. Defaults to 1.
    
Returns:
    float: Bayesian AIC 
'''
```
```Python
DOF(ndat, nparam, priorsigma=None) -> int:
'''
Compute the number of degrees of freedom. Depends on whether you use priors. Any input priors are taken as
initial guesses for the fit algorithm. If you would like parameters in the prior array to be treated as a 
starting guess only, and not as a Bayesian prior, set its corresponding error to np.inf. Hence when there
are priors, the number of degrees of freedom equals the number of ydata, less the number of finite prior errors.

Args:
    ndat (int): number of data 
    nparam (int): number of model parameters 
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

Returns:
    int: number of degrees of freedom 
'''
```
```Python
KSTest_1side(data, cdf) -> float:
'''
1-sided Kolmogorov test. Gives back the likelihood that the observed difference between
data and cdf are at least as extreme as suggested by the Kolmogorov statistic.

Args:
    data (np.ndarray)
    cdf (function)

Returns:
    float: 1-p 
'''
```
```Python
KSTest_2side(data1, data2) -> float:
'''
2-sided Kolmogorov test. Gives back the likelihood that the observed difference between
data1 and data2 are at least as extreme as suggested by the Kolmogorov statistic.

Args:
    data1 (np.ndarray)
    data2 (np.ndarray)

Returns:
    float: 1-p 
'''
```
```Python
biased_sample_variance(data, err) -> float:
'''
Compute the biased weighted sample variance, i.e. the biased variance of an 
individual measurement and not the variance of the mean.

Args:
    data (array-like)
    err (array-like)

Returns:
    float: sample variance 
'''
```
```Python
binSeries(data, nbins) -> numpy.ndarray:
'''
Take a time series and bin it. Bin 0 is the average over the first binsize elements,
bin 1 the average over the next binsize elements, and so on.

Args:
    data (array-like)
    nbins (int)

Returns:
    np.ndarray: Binned data
'''
```
```Python
checkDomain(domain):
'''
Some methods require that you do something over an interval, which we refer to in this module as
a 'domain'. This checks the domain makes sense.

Args:
    domain (tuple)
'''
```
```Python
checkPrior(prior, priorsigma):
'''
Make sure prior and priorsigma status are compatible. 
'''
```
```Python
checkProb(p)
```
```Python
checkTS(ts):
'''
Some methods require 1-d time series. This checks that the type, dimensionality,
and length are appropriate.

Args:
    ts (array-like): time series 
'''
```
```Python
chisquare(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
'''
Calculate chi^2, see e.g. eq. (8.28) of Sivia and Skilling or eq. (A1) of
10.1103/PhysRevD.90.054506. We assume priors are not correlated with data.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): model parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
    
    Returns:
        float: chi^2
'''
```
```Python
confidence_ellipse(x, y, CI=None, **params):
'''
Plot a confidence ellipse according to the data x, y. The confidence is only meaningful 
assuming the x and y are Gaussian distributed. By default, draws an ellipse that captures
roughly 39% of the data.

Args:
    x (array-like)
    y (array-like)
    ax (matplotlib ax object)
    color (str, optional): Color of the ellipse edge. Defaults to 'r'.
    C (float, optional): Desired confidence. Defaults to ~0.39.

Returns:
    float, float: semi-major and semi-minor lengths of drawn ellipse 
'''
```
```Python
countParams(func, params) -> int:
'''
Count number of model parameters. For a typical function without priors,
we count the length of the params array. Otherwise we assume it's a spline.

Args:
    func (func)
    params (array-like): model parameters.
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

Returns:
    int: number of parameters. 
'''
```
```Python
countPriors(priorsigma=None) -> int:
'''
The number of priors is the number of finite prior error bars.

Args:
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

Returns:
    int: Number of priors 
'''
```
```Python
cov_to_cor(cov) -> numpy.ndarray:
'''
Normalize covariance matrix to create correlation matrix.

Args:
    cov (np.ndarray)

Returns:
    np.ndarray: correlation matrix 
'''
```
```Python
covariance(x, y) -> float:
'''
Unbiased estimator of the covariance between the time series x and y.

Args:
    x (array-like)
    y (array-like)

Returns:
    float: cov 
'''
```
```Python
dev_by_dist(data, axis=0, return_both_q=False, percentile=68):
'''
Calculate the distance between the median and 68% quantiles. Returns the larger of the two 
distances. This method is used sometimes to estimate error, for example in the bootstrap. 
'''
```
```Python
empiricalCDF(data):
'''
Create the x and y coordinates needed to plot the empirical CDF
of a 1-d set of data.

Args:
    data (array-like): measurements 

Returns:
    func: CDF 
'''
```
```Python
error_prop(func, means, errors, grad=None, args=()):
'''
Use error propagation to propagate some errors through function func. The function should have the form
    func( data ), 
where data is a 1-d array of input variables. 

Args:
    func (func)
    means (np.ndarray)
    errors (np.ndarray)
    grad (func, optional): Gradient function. Defaults to None.
    args (tuple, optional): Arguments of func. Defaults to ().

Returns:
    np.ndarray, np.ndarray: f, f_err 
'''
```
```Python
error_prop_func(x, func, params, params_err, grad=None, args=()):
'''
Propagate error in f(x;params,params_err). This needs its own special treatment, since
the error propagation method on its own only propagates params_err to f(params,params_err).

Args:
    x (array-like)
    func (func)
    params (array-like): Model parameters. 
    params_err (array-like): Error in model parameters. 
    grad (func, optional): Gradient function. Defaults to None.
    args (tuple, optional): Arguments of func. Defaults to ().
'''
```
```Python
expandArgs(func, x, params=(), args=()):
'''
In general we distinguish between parameters and arguments. Parameters should be passed
together as a collection, e.g. as a tuple, list, or np.array. Other function arguments can
be passed how you like and will be expanded here.

Args:
    func (func)
    x (array-like)
    params (tuple, optional): Model parameters. Defaults to ().
    args (tuple, optional): Function arguments. Defaults to ().

Returns:
    func(x,params,args) 
'''
```
```Python
forcePositiveSemidefinite(mat):
'''Doctors a noisy correlation matrix mat to be positive semidefinite if it isn't already. 
Uses algorithm of Rebonato and Jaeckel, DOI: 10.2139/ssrn.1969689

Args:
    mat (np.ndarray)

Returns:
    np.ndarray: positive semidefinite matrix 
'''
```
```Python
gaudif(x1, e1, x2, e2) -> float:
'''
Likelihood that difference between outcomes x1 and x2 is due to chance, assuming x1 and x2 are
both drawn from a normal distribution with the same mean. A rule of thumb is that this is more
appropriate when one estimated x1 and x2 using ~30 or more measurements.

Args:
    x1 (float): mean 1 
    e1 (float): standard error 1
    x2 (float): mean 2
    e2 (float): standard error 2

Returns:
    float: p-value 
'''
```
```Python
getModelWeights(IC) -> numpy.ndarray:
'''
Convert information criteria IC to normalized probability weights.

Args:
    IC (array-like): Array of information criteria 

Returns:
    np.array: Probability weights 
'''
```
```Python
goodnessOfFit(dof, chi2) -> float:
'''
The q-value or goodness of fit.

Args:
    dof (int): number of degrees of freedom 
    chi2 (float): the chi2 

Returns:
    float: Q 
'''
```
```Python
logGBF(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
'''
log P(data|model). This quantity is useful for comparing fits of the same data to different models that
have different priors and/or fit functions. The model with the largest logGBF is the one preferred by the data.
Differences in logGBF smaller than 1 are not very significant. Gaussian statistics are assumed.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): model parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
    
    Returns:
        float: log( Gaussian Bayes factor )
'''
```
```Python
meanArgWrapper(func, used_data, args)
```
```Python
midpointMeanError(lo, hi):
'''
Take mean+err and mean-err and use those to compute mean and err. Assume symmetric
error bar. This is useful e.g. when using WebPlotDigitizer.

Args:
    lo (float): mean-err 
    hi (float): mean+err 

Returns:
    mean, err 
'''
```
```Python
modelAverage(data, err, IC, return_syst=False):
'''
Given some fit results, corresponding error, and information criteria, compute
a weighted model average.

Args:
    data (array-like): Fit results 
    err (array-like): Result errors 
    IC (array-like): Information criteria 

Returns:
    tuple: Model average and error (optionally systematic error)
'''
```
```Python
pearson(x, y) -> float:
'''
Get the Pearson correlation coefficient between the time series x and y.

Args:
    x (array-like)
    y (array-like)

Returns:
    float: R 
'''
```
```Python
plot_correlation(mat, ax=<module 'matplotlib.pyplot' from '/home/dclarke/.local/lib/python3.14/site-packages/matplotlib/pyplot.py'>):
'''Plot correlation matrix as a heatmap.

Args:
    mat (np.ndarray): correlation matrix
    ax (matplotlib ax object): Defaults to plt.
'''
```
```Python
plot_func(func, domain, params=(), args=(), func_err=None, params_err=(), grad=None, swapXY=False, npoints=1000, **kwargs):
'''
Plot a function along with its error bands.

Args:
    func (func)
    domain (tuple): Domain of function. 
    params (tuple, optional): Model parameters. Defaults to ().
    params_err (tuple, optional): Error in model parameters. Defaults to ().
    args (tuple, optional): Optional function arguments. Defaults to ().
    func_err (func, optional): Explicit error function. Defaults to None.
    grad (func, optional): Explicit function gradient to compute error. Defaults to None.
    swapXY (bool, optional): Swap X and Y variables in plot. Defaults to False.
    npoints (int, optional): Number of points to use for plotting. Defaults to 1000.
'''
```
```Python
save_func(func, domain, params=(), args=(), func_err=None, params_err=(), grad=None, npoints=1000, header=None, filename='func.d', **kwargs):
'''
Save a function along with its error bands.

Args:
    func (func)
    domain (tuple): Domain of function. 
    params (tuple, optional): Model parameters. Defaults to ().
    params_err (tuple, optional): Error in model parameters. Defaults to ().
    args (tuple, optional): Optional function arguments. Defaults to ().
    func_err (func, optional): Explicit error function. Defaults to None.
    grad (func, optional): Explicit function gradient to compute error. Defaults to None.
    swapXY (bool, optional): Swap X and Y variables in plot. Defaults to False.
    npoints (int, optional): Number of points to use for plotting. Defaults to 1000.
'''
```
```Python
std_dev(data, axis=0):
'''
Calculate unbiased (ddof = 1) estimator for the standard deviation. 
 
    The default behavior of numpy is to flatten the data, flagged by axis=None. This is
    something that is never needed in our context. Changing the default to axis=0 means
    applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
    '''
```
```Python
std_err(data, axis=0):
'''
Standard deviation of the sample mean according to the CLT. 
 
    The default behavior of numpy is to flatten the data, flagged by axis=None. This is
    something that is never needed in our context. Changing the default to axis=0 means
    applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
    '''
```
```Python
std_mean(data, axis=0):
'''
Compute the mean. 
 
    The default behavior of numpy is to flatten the data, flagged by axis=None. This is
    something that is never needed in our context. Changing the default to axis=0 means
    applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
    '''
```
```Python
std_median(data, axis=0):
'''
Compute the median. 
 
    The default behavior of numpy is to flatten the data, flagged by axis=None. This is
    something that is never needed in our context. Changing the default to axis=0 means
    applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
    '''
```
```Python
std_var(data, axis=0):
'''
Calculate unbiased (ddof = 1) estimator for the variance. 
 
    The default behavior of numpy is to flatten the data, flagged by axis=None. This is
    something that is never needed in our context. Changing the default to axis=0 means
    applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
    '''
```
```Python
studif(x1, e1, ndat1, x2, e2, ndat2) -> float:
'''
Likelihood that difference between outcomes x1 and x2 is due to chance, assuming x1 and x2 are
both drawn from a normal distribution with the same mean. A rule of thumb is that this is more
appropriate when one estimated x1 and x2 using ~30 or fewer measurements. Of course, you can
always compare this with gaudif to get a better idea.

Args:
    x1  (float): mean 1 
    e1  (float): standard error 1
    ndat1 (int): number of measurements used to compute x1
    x2  (float): mean 2
    e2  (float): standard error 2
    ndat2 (int): number of measurements used to compute x2

Returns:
    float: p-value 
'''
```
```Python
symmetrizeError(lo, hi, central, method='conservative'):
'''
Take unsymmetric errors and estimate symmetric errors from them.

Args:
    lo (float): mean-err1 
    hi (float): mean+err2
    central (float): mean 
    method (str, optional): Approach to estimate error. Defaults to 'conservative'.

Returns:
    mean and symmetric error
'''
```
```Python
unbiased_mean_variance(data, err) -> float:
'''
Compute the unbiased variance of a weighted mean. Do not use this function if your weights are frequency
weights. This is more like a systematic error. The absolute size of the weights does not matter. The error is
constructed using the deviations of the individual data points. 
'''
```
```Python
unbiased_sample_variance(data, err) -> float:
'''
Compute the unbiased weighted sample variance, i.e. the unbiased variance of an individual measurement and not
the variance of the mean. Do not use this function if your weights are frequency weights. 
'''
```
```Python
weighted_mean(data, err) -> float:
'''
Compute the weighted mean. Here the weights are Gaussian error bars.
See e.g. https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html.

Args:
    data (array-like)
    errcov (array-like): error if 1d or cov if 2d. 

Returns:
    float: weighted mean 
'''
```
```Python
weighted_variance(err) -> float:
'''
Get variance of above weighted mean, when the weights are statistical errors. 

Args:
    err (array-like)

Returns:
    float: weighted variance 
'''
```
