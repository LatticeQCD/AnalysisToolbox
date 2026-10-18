latqcdtools.statistics.statistics
=============

```python
AIC(xdata, ydata, cov, func, args=(), params=()) -> float:
"""
The Akaike information criterion (AIC) is a measure of how well a fit performs. It builds on the likelihood
function by including a penalty for each d.o.f. This is useful in a context where you have multiple models to
choose from,and hence different numbers of d.o.f. possible. It's also useful when you are worried about
overfitting. The preferred model minimizes the AIC. With k the number of parameters,
    AIC = 2k - 2 log L = chi^2 + 2k + log det C + N log(2 pi).
The chi^2 of a fit is lower than the chi^2 the same fit would have on fresh data, since the fit adapts to
the noise. The 2k corrects this bias. For a model linear in its parameters and a known C, the bias is
exactly 2k for any N. For a smoothing spline, k is the effective number tr(H). The last two terms are the
same for every model of the same data, so they drop out of comparisons.

The AIC does not know about priors: it uses the data-only chi^2 and counts every parameter fully. For a
fit with priors it is therefore not the right criterion. Use logGBF or BAIC instead.

Args:
    xdata (array-like)
    ydata (array-like)
    cov (array-like): covariance matrix 
    func (func)
    args (tuple, optional): arguments to func. Defaults to ().
    params (tuple, optional): model parameters. Defaults to ().

Returns:
    float: AIC
"""
```
```python
AICc(xdata, ydata, cov, func, args=(), params=()) -> float:
"""
Corrected AIC (AICc) of Sugiura, Commun. Stat. Theory Methods 7, 13 (1978) and Hurvich & Tsai, Biometrika 76,
297 (1989); see also Burnham & Anderson, "Model Selection and Multimodel Inference" (2nd ed., Springer, 2002),
Eq. (2.4). With k the number of parameters and N the number of data,
    AICc = AIC + 2k(k+1)/(N-k-1).
The extra term corrects the small-sample bias that comes from estimating the noise variance from the
residuals. Here the covariance C is given rather than estimated from the residuals, and then the AIC
is already unbiased for a model linear in its parameters, so the AICc over-penalizes parameters at small N.
Prefer AIC. Like the AIC, the AICc does not know about priors, so it is not the right criterion for a fit
with priors. For N <= k+1 the correction is undefined, and we return +inf, so that in a model average
(weight ~ exp(-AICc/2)) such a model gets zero weight.

Args:
    xdata (array-like)
    ydata (array-like)
    cov (array-like): covariance matrix 
    func (func)
    args (tuple, optional): arguments to func. Defaults to ().
    params (tuple, optional): model parameters. Defaults to ().

Returns:
    float: corrected AIC
"""
```
```python
BAIC(xdata, ydata, cov, func, args=(), params=(), Ncut=0, modelPrior=1) -> float:
"""
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
"""
```
```python
DOF(ndat, nparam, priorsigma=None) -> int:
"""
Compute the number of degrees of freedom. Depends on whether you use priors. Any input priors are taken as
initial guesses for the fit algorithm. If you would like parameters in the prior array to be treated as a 
starting guess only, and not as a Bayesian prior, set its corresponding error to np.inf. Hence when there
are priors, the number of degrees of freedom equals the number of ydata, less the number of finite prior errors.

Args:
    ndat (int): number of data 
    nparam (int): number of model parameters. May be an effective, non-integer number, e.g. for a 
      smoothing spline.
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

Returns:
    int: number of degrees of freedom 
"""
```
```python
KSTest_1side(data, cdf) -> float:
"""
1-sided Kolmogorov test. Gives back the likelihood that the observed difference between
data and cdf are at least as extreme as suggested by the Kolmogorov statistic.

Args:
    data (np.ndarray)
    cdf (function)

Returns:
    float: 1-p 
"""
```
```python
KSTest_2side(data1, data2) -> float:
"""
2-sided Kolmogorov test. Gives back the likelihood that the observed difference between
data1 and data2 are at least as extreme as suggested by the Kolmogorov statistic.

Args:
    data1 (np.ndarray)
    data2 (np.ndarray)

Returns:
    float: 1-p 
"""
```
```python
_inverseVarianceWeights(err) -> numpy.ndarray:
"""
Weights 1/err^2 for the weighted-mean methods below. 
"""
```
```python
biased_sample_variance(data, err) -> float:
"""
Compute the biased weighted sample variance, i.e. the biased variance of an 
individual measurement and not the variance of the mean.

Args:
    data (array-like)
    err (array-like)

Returns:
    float: sample variance 
"""
```
```python
binSeries(data, nbins) -> numpy.ndarray:
"""
Take a time series and bin it. Bin 0 is the average over the first binsize elements,
bin 1 the average over the next binsize elements, and so on.

Args:
    data (array-like)
    nbins (int)

Returns:
    np.ndarray: Binned data
"""
```
```python
checkDomain(domain):
"""
Some methods require that you do something over an interval, which we refer to in this module as
a 'domain'. This checks the domain makes sense.

Args:
    domain (tuple)
"""
```
```python
checkPrior(prior, priorsigma):
"""
Make sure prior and priorsigma status are compatible. 
"""
```
```python
checkProb(p)
```
```python
checkTS(ts):
"""
Some methods require 1-d time series. This checks that the type, dimensionality,
and length are appropriate.

Args:
    ts (array-like): time series 
"""
```
```python
chisquare(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
"""
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
"""
```
```python
confidence_ellipse(x, y, CI=None, **params):
"""
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
"""
```
```python
countParams(func, params) -> int:
"""
Count number of model parameters. For a typical function without priors,
we count the length of the params array. Otherwise we assume it's a spline.

Args:
    func (func)
    params (array-like): model parameters.
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

Returns:
    int: number of parameters. 
"""
```
```python
countPriors(priorsigma=None) -> int:
"""
The number of priors is the number of finite prior error bars.

Args:
    priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

Returns:
    int: Number of priors 
"""
```
```python
cov_to_cor(cov) -> numpy.ndarray:
"""
Normalize covariance matrix to create correlation matrix.

Args:
    cov (np.ndarray)

Returns:
    np.ndarray: correlation matrix 
"""
```
```python
covariance(x, y) -> float:
"""
Unbiased estimator of the covariance between the time series x and y.

Args:
    x (array-like)
    y (array-like)

Returns:
    float: cov 
"""
```
```python
dev_by_dist(data, axis=0, return_both_q=False, percentile=68):
"""
Calculate the distance between the median and the quantiles enclosing the central percentile of the
data. This method is used sometimes to estimate error, for example in the bootstrap. The quantiles
are not interpolated; they are rounded outward to the nearest datum. If many data share the median
value, both distances can be zero even though the data have a spread.

Args:
    data (array-like)
    axis (int, optional): Axis along which to compute. Defaults to 0.
    return_both_q (bool, optional): Return both distances instead of the larger one? Defaults to False.
    percentile (float, optional): Central percentile. Defaults to 68.

Returns:
    distance to the lower and upper quantiles, or the larger of the two
"""
```
```python
empiricalCDF(data):
"""
Create the x and y coordinates needed to plot the empirical CDF
of a 1-d set of data.

Args:
    data (array-like): measurements 

Returns:
    func: CDF 
"""
```
```python
error_prop(func, means, errors, grad=None, args=()):
"""
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
"""
```
```python
error_prop_func(x, func, params, params_err, grad=None, args=()):
"""
Propagate error in f(x;params,params_err). This needs its own special treatment, since
the error propagation method on its own only propagates params_err to f(params,params_err).

Args:
    x (array-like)
    func (func)
    params (array-like): Model parameters. 
    params_err (array-like): Covariance matrix of model parameters
    grad (func, optional): Gradient function. Defaults to None.
    args (tuple, optional): Arguments of func. Defaults to ().
"""
```
```python
estimateCovariance(data_std_dev) -> numpy.ndarray:
"""
Covariance matrix of some data, under the assumption that each datum X_i is maximally correlated
with every other datum X_j, i.e. cov_ij = std_i std_j. This matrix has rank 1, so it is not invertible.

Args:
    data_std_dev (array-like): One-dimensional array of standard deviations.

Returns:
    np.ndarray: covariance matrix
"""
```
```python
expandArgs(func, x, params=(), args=()):
"""
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
"""
```
```python
forcePositiveSemidefinite(mat):
"""Doctors a noisy correlation matrix mat to be positive semidefinite if it isn't already. 
Uses algorithm of Rebonato and Jaeckel, DOI: 10.2139/ssrn.1969689

Args:
    mat (np.ndarray)

Returns:
    np.ndarray: positive semidefinite matrix 
"""
```
```python
gaudif(x1, e1, x2, e2) -> float:
"""
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
"""
```
```python
getModelWeights(IC) -> numpy.ndarray:
"""
Convert information criteria IC to normalized probability weights.

Args:
    IC (array-like): Array of information criteria 

Returns:
    np.array: Probability weights 
"""
```
```python
goodnessOfFit(dof, chi2) -> float:
"""
The q-value or goodness of fit.

Args:
    dof (int): number of degrees of freedom 
    chi2 (float): the chi2 

Returns:
    float: Q 
"""
```
```python
logGBF(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None, jac=None) -> float:
"""
log P(data|model). This quantity is useful for comparing fits of the same data to different models that
have different priors and/or fit functions. The model with the largest logGBF is the one preferred by the data.
Differences in logGBF smaller than 1 are not very significant. Gaussian statistics are assumed.

The evidence is the integral over parameters of likelihood times prior. Expanding chi^2 to second order
about its minimum (the Laplace approximation, exact for a linear model) gives
    logGBF = 1/2 [ log det Sigma_post - log det C - log det Sigma_prior - chi^2 - dof log(2 pi) ],
with C the data covariance, chi^2 the augmented chi^2 at params, and
    Sigma_post^{-1} = J^T C^{-1} J + Sigma_prior^{-1},   J_ij = d f(x_i) / d p_j.
The J^T C^{-1} J term drops second derivatives of f, as in lsqfit. A parameter with infinite priorsigma
gets a flat prior of unit density: it drops out of Sigma_prior and the prior count in dof. With no priors
at all, this is the flat-prior evidence. Each factor of 2 pi comes from one data point, one finite prior,
or (with the opposite sign) one integrated parameter, so the 2 pi factors combine to dof log(2 pi). When
every parameter has a finite prior, dof equals the number of data, and this matches logGBF of lsqfit.
If jac is not given, J is computed with numerical derivatives. A Fitter passes the J it used for its
parameter covariance, so that Sigma_post and pcov agree.

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
"""
```
```python
logLikelihood(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
"""
Log of the Gaussian likelihood P(data|params) at params, which is usually the best fit:
    log L = -1/2 [ chi^2 + log det C + N log(2 pi) ],
with C the data covariance and N the number of data. With priors, each finite prior is treated as one
more Gaussian datum, so that L is the likelihood times the prior, P(data|params) P(params):
    log L = -1/2 [ chi^2 + log det C + log det Sigma_prior + (N + nprior) log(2 pi) ],
with chi^2 the augmented chi^2 and Sigma_prior = diag(priorsigma^2) over the nprior finite priors.

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
    float: log of the likelihood
"""
```
```python
meanArgWrapper(func, used_data, args)
```
```python
midpointMeanError(lo, hi):
"""
Take mean+err and mean-err and use those to compute mean and err. Assume symmetric
error bar. This is useful e.g. when using WebPlotDigitizer.

Args:
    lo (float): mean-err 
    hi (float): mean+err 

Returns:
    mean, err 
"""
```
```python
modelAverage(data, err, IC, return_syst=False):
"""
Given some fit results, corresponding error, and information criteria, compute
a weighted model average.

Args:
    data (array-like): Fit results 
    err (array-like): Result errors 
    IC (array-like): Information criteria 

Returns:
    tuple: Model average and error (optionally systematic error)
"""
```
```python
pearson(x, y) -> float:
"""
Get the Pearson correlation coefficient between the time series x and y.

Args:
    x (array-like)
    y (array-like)

Returns:
    float: R 
"""
```
```python
plot_correlation(mat, ax=<module 'matplotlib.pyplot' from '/Users/dclarke/GitHub/AnalysisToolbox/venv/AnalysisToolbox/lib/python3.9/site-packages/matplotlib/pyplot.py'>, vmin=-1, vmax=1):
"""Plot correlation matrix as a heatmap.

Args:
    mat (np.ndarray): correlation matrix
    ax (matplotlib ax object): Defaults to plt.
"""
```
```python
plot_func(func, domain, params=(), args=(), func_err=None, params_err=(), grad=None, swapXY=False, npoints=1000, **kwargs):
"""
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
"""
```
```python
save_func(func, domain, params=(), args=(), func_err=None, params_err=(), grad=None, npoints=1000, header=None, filename='func.d', **kwargs):
"""
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
"""
```
```python
std_dev(data, axis=0):
"""
Calculate unbiased (ddof = 1) estimator for the standard deviation. 
 
The default behavior of numpy is to flatten the data, flagged by axis=None. This is
something that is never needed in our context. Changing the default to axis=0 means
applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
"""
```
```python
std_err(data, axis=0):
"""
Standard deviation of the sample mean according to the CLT. 
 
The default behavior of numpy is to flatten the data, flagged by axis=None. This is
something that is never needed in our context. Changing the default to axis=0 means
applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
"""
```
```python
std_mean(data, axis=0):
"""
Compute the mean. 
 
The default behavior of numpy is to flatten the data, flagged by axis=None. This is
something that is never needed in our context. Changing the default to axis=0 means
applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
"""
```
```python
std_median(data, axis=0):
"""
Compute the median. 
 
The default behavior of numpy is to flatten the data, flagged by axis=None. This is
something that is never needed in our context. Changing the default to axis=0 means
applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
"""
```
```python
std_var(data, axis=0):
"""
Calculate unbiased (ddof = 1) estimator for the variance. 
 
The default behavior of numpy is to flatten the data, flagged by axis=None. This is
something that is never needed in our context. Changing the default to axis=0 means
applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
"""
```
```python
studif(x1, e1, ndat1, x2, e2, ndat2) -> float:
"""
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
"""
```
```python
symmetrizeError(lo, hi, central, method='conservative'):
"""
Take unsymmetric errors and estimate symmetric errors from them.

Args:
    lo (float): mean-err1 
    hi (float): mean+err2
    central (float): mean 
    method (str, optional): Approach to estimate error. Defaults to 'conservative'.

Returns:
    mean and symmetric error
"""
```
```python
unbiased_mean_variance(data, err) -> float:
"""
Compute the unbiased variance of a weighted mean. Do not use this function if your weights are frequency
weights. This is more like a systematic error. The absolute size of the weights does not matter. The error is
constructed using the deviations of the individual data points. This assumes all data share one unknown
variance; err only sets the relative weights 1/err^2 (reliability weights). 
"""
```
```python
unbiased_sample_variance(data, err) -> float:
"""
Compute the unbiased weighted sample variance, i.e. the unbiased variance of an individual measurement and not
the variance of the mean. Do not use this function if your weights are frequency weights. This assumes all
data share one unknown variance; err only sets the relative weights 1/err^2 (reliability weights). 
"""
```
```python
weighted_mean(data, err) -> float:
"""
Compute the weighted mean. Here the weights are Gaussian error bars.
See e.g. https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html.

Args:
    data (array-like)
    errcov (array-like): error if 1d or cov if 2d. 

Returns:
    float: weighted mean 
"""
```
```python
weighted_variance(err) -> float:
"""
Get variance of above weighted mean, when the weights are statistical errors. 

Args:
    err (array-like)

Returns:
    float: weighted variance 
"""
```
```python
whiteningMatrix(cov, svdcut=None) -> numpy.ndarray:
"""
Whitening matrix W of a covariance matrix C, i.e. a matrix satisfying W^T W = C^-1.

Let r = f(p) - y be the residuals of a fit. Their whitened version z = W r has two properties:
    1. chi^2 = r^T C^-1 r = z^T z, i.e. chi^2 is a plain sum of squares;
    2. Cov(z) = W C W^T = 1, i.e. the components of z are uncorrelated with unit variance. 
For uncorrelated data, W = diag(1/sigma), so z_i = r_i/sigma_i. For correlated data, W also rotates the
residuals into the eigenbasis of the correlation matrix. Whitening turns a correlated fit into an
ordinary least-squares problem, which lets a least-squares solver use the Jacobian dz/dp = W df/dp to take
Gauss-Newton steps. These converge quickly and precisely near the minimum. 

Many matrices satisfy W^T W = C^-1. Following gvar/lsqfit, we diagonalize the correlation matrix
cor = V Lambda V^T and use W = Lambda^(-1/2) V^T diag(1/sigma); see CovarianceDecomposition. Small
eigenvalues of cor get inverted here, so an ill-conditioned cor needs an SVD cut. 

Args:
    cov (np.ndarray): covariance matrix 
    svdcut (float, optional): SVD cut applied to the correlation matrix. Defaults to None, meaning no cut.

Returns:
    np.ndarray: whitening matrix W 
"""
```
```python
class CovarianceDecomposition(cov, svdcut=None):
"""
Decomposition of a covariance matrix C, following gvar/lsqfit. Write C = D cor D, with D = diag(sigma) and
cor the correlation matrix, and diagonalize cor = V Lambda V^T. Everything a fit needs follows from this one
decomposition:
    whitening matrix   W = Lambda^(-1/2) V^T D^(-1),
    inverse            C^(-1) = W^T W,
    log determinant    log det C = 2 sum_i log sigma_i + sum_i log lambda_i.
The SVD cut is applied to the eigenvalues of cor, so that it does not depend on the units of the data. 
For a positive semidefinite matrix this is the same as regulate(cor). 

Args:
    cov (np.ndarray): covariance matrix 
    svdcut (float, optional): Raise eigenvalues of cor below svdcut*max to that floor. Defaults to None,
      meaning no cut.
"""
```
