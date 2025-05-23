#
# statistics.py
#
# D. Clarke, H Sandmeyer 
#
# A collection of basic methods for statistical analysis. The following methods have been adapted from software of
# Bernd Berg, Markov Chain Monte Carlo and their Statistical Analysis, World Scientific, 2004, ISBN=978-981-3106-37-6:
#     gaudif, studif
#


import numpy as np
import scipy as sp
from latqcdtools.math.num_deriv import diff_jac 
from latqcdtools.math.math import logDet, normalize, invert
from latqcdtools.base.plotting import fill_param_dict, plot_fill, plot_lines, FOREGROUND, plt
from latqcdtools.base.utilities import isHigherDimensional, toNumpy, appendToDocstring, unvector,\
    createFilePath
from latqcdtools.base.cleanData import clipRange
from latqcdtools.base.check import checkType, checkEqualLengths
from latqcdtools.base.readWrite import writeTable
import latqcdtools.base.logger as logger
import matplotlib


NUMPYCOMMENT = """ 
    The default behavior of numpy is to flatten the data, flagged by axis=None. This is
    something that is never needed in our context. Changing the default to axis=0 means
    applying this function to an np.ndarray of shape (N,M) yields an array of shape (M,). 
    """

ICARGCOMMENTS = """
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): model parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
    """


def meanArgWrapper(func,used_data,args):
    if isinstance(args, dict):
        return func(used_data, **args)
    else:
        return func(used_data, *args)


@appendToDocstring(NUMPYCOMMENT)
def std_median(data, axis = 0):
    """ 
    Compute the median. 
    """
    return np.median(data, axis)


@appendToDocstring(NUMPYCOMMENT)
def std_mean(data, axis = 0):
    """ 
    Compute the mean. 
    """
    return np.mean(data, axis)


@appendToDocstring(NUMPYCOMMENT)
def std_var(data, axis = 0):
    """ 
    Calculate unbiased (ddof = 1) estimator for the variance. 
    """
    data = np.asarray(data)
    return np.var(data, axis = axis, ddof = 1)


@appendToDocstring(NUMPYCOMMENT)
def std_dev(data, axis = 0):
    """ 
    Calculate unbiased (ddof = 1) estimator for the standard deviation. 
    """
    data = np.asarray(data)
    return np.std(data, axis = axis, ddof = 1)


@appendToDocstring(NUMPYCOMMENT)
def std_err(data, axis = 0):
    """ 
    Standard deviation of the sample mean according to the CLT. 
    """
    data = np.asarray(data)
    return std_dev(data, axis) / np.sqrt(data.shape[axis])


def expandArgs(func, x, params=(), args=()):
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
    if len(params)==0:
        # Python treats func(x,*()) as func(x)
        return func(x, *args)
    else:
        return func(x, params, *args)


def checkDomain(domain):
    """ 
    Some methods require that you do something over an interval, which we refer to in this module as
    a 'domain'. This checks the domain makes sense.

    Args:
        domain (tuple)
    """
    checkType(tuple,domain=domain)
    if len(domain) != 2:
        logger.TBRaise('A domain is a tuple of the form (xmin,xmax) specifying the interval [xmin,xmax].',frame=3)
    if domain[0]>=domain[1]:
        logger.TBRaise('Must have domain[1]>domain[0].',frame=3)


def checkPrior(prior,priorsigma):
    """ 
    Make sure prior and priorsigma status are compatible. 
    """
    if prior is None and priorsigma is not None:
        logger.TBRaise('prior = None, priorsigma != None',frame=3)
    if priorsigma is None and prior is not None:
        logger.TBRaise('prior != None, priorsigma = None',frame=3)


def checkTS(ts):
    """ 
    Some methods require 1-d time series. This checks that the type, dimensionality,
    and length are appropriate.

    Args:
        ts (array-like): time series 
    """
    checkType(np.ndarray,ts=ts)
    if isHigherDimensional(ts):
        logger.TBRaise('Expected 1-d time series.',frame=3)
    if len(ts) < 2:
        logger.TBRaise('Time series needs at least two measurements.',frame=3)


def checkProb(p):
    checkType('real',p=p)
    if not 0 <= p <= 1:
        logger.TBRaise('Probabilities must be between 0 and 1.',frame=3)


def countParams(func,params) -> int:
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
    nparam = len(params)
    if nparam == 0:
        try:
            # customSpline
            nparam = len(func.get_coeffs())
            # CubicSpline
        except AttributeError:
            nparam = len(func.c)
    return nparam


def countPriors(priorsigma=None) -> int:
    """ 
    The number of priors is the number of finite prior error bars.

    Args:
        priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.

    Returns:
        int: Number of priors 
    """
    nprior = 0
    if priorsigma is not None:
        nprior = np.count_nonzero(priorsigma!=np.inf)
    return nprior


def DOF(ndat,nparam,priorsigma=None) -> int:
    """  
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
    """
    checkType(int,ndat=ndat)
    checkType(int,nparam=nparam)
    nprior = countPriors(priorsigma) 
    dof = ndat + nprior - nparam 
    logger.debug('dof =',dof,'ndat =',ndat,'nparam =',nparam,'nprior =',nprior)
    if dof < 0:
        logger.TBRaise('ndat =',ndat,'nparam =',nparam,'nprior =',nprior)
    return dof


@appendToDocstring(args=ICARGCOMMENTS,returns='\n        float: chi^2\n')
def chisquare(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
    """ 
    Calculate chi^2, see e.g. eq. (8.28) of Sivia and Skilling or eq. (A1) of
    10.1103/PhysRevD.90.054506. We assume priors are not correlated with data.
    """
    checkPrior(prior,priorsigma)
    y    = expandArgs(func,xdata,params,args)
    diff = ydata-y
    res  = diff @ invert(cov) @ diff
    if prior is not None:
        res += np.sum((np.array(params) - prior)**2 / priorsigma**2)
    return res


@appendToDocstring(args=ICARGCOMMENTS,returns='\n        float: log( Gaussian Bayes factor )\n' ) 
def logGBF(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
    """ 
    log P(data|model). This quantity is useful for comparing fits of the same data to different models that
    have different priors and/or fit functions. The model with the largest logGBF is the one preferred by the data.
    Differences in logGBF smaller than 1 are not very significant. Gaussian statistics are assumed.
    """
    chi2   = chisquare(xdata, ydata, cov, func, args, params, prior, priorsigma)
    nparam = countParams(func,params)
    dof    = DOF(len(ydata),nparam,priorsigma)
    if prior is None:
        return 0.5*( - logDet(cov) - chi2 - dof*np.log(2*np.pi) )
    else:
        return 0.5*( - logDet(cov) - chi2 - dof*np.log(2*np.pi) + logDet(np.diag(clipRange(priorsigma)**2)) )


@appendToDocstring(args=ICARGCOMMENTS,returns='\n        float: AIC\n')
def AIC(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
    """ 
    The Akaike information criterion (AIC) is a measure of how well a fit performs. It builds on the likelihood
    function by including a penalty for each d.o.f. This is useful in a context where you have multiple models to
    choose from,and hence different numbers of d.o.f. possible. It's also useful when you are worried about
    overfitting. The preferred model minimizes the AIC.
    """
    nparam     = countParams(func,params)
    likelihood = logGBF(xdata, ydata, cov, func, args, params, prior, priorsigma)
    return 2*nparam - 2*likelihood


@appendToDocstring(args=ICARGCOMMENTS,returns='\n        float: corrected AIC\n')
def AICc(xdata, ydata, cov, func, args=(), params=(), prior=None, priorsigma=None) -> float:
    """ 
    Corrected AIC (AICc). When the sample size is smaller, it increases the chance AIC will select a model with too
    many parameters. The AICc tries to further correct for this. In the limit that the number of data points goes to
    infinity, one recovers the AIC.
    """
    nparam = countParams(func,params)
    nprior = countPriors(priorsigma) 
    ndat   = len(ydata) + nprior
    aic    = AIC(xdata, ydata, cov, func, args, params, prior, priorsigma)
    return aic + 2*(nparam**2+nparam)/(ndat-nparam+1)


def BAIC(xdata, ydata, cov, func, args=(), params=(), Ncut=0, modelPrior=1) -> float:
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
    nparam   = countParams(func,params)
    chi2data = chisquare(xdata, ydata, cov, func, args, params, None, None)
    logger.debug('chi_data^2 =',chi2data)
    return -2*np.log(modelPrior) + chi2data + 2*nparam + 2*Ncut 


def pearson(x,y) -> float:
    """ 
    Get the Pearson correlation coefficient between the time series x and y.

    Args:
        x (array-like)
        y (array-like)

    Returns:
        float: R 
    """
    checkTS(x)
    checkTS(y)
    checkEqualLengths(x,y)
    return np.corrcoef(x,y)[0,1]


def covariance(x,y) -> float:
    """ 
    Unbiased estimator of the covariance between the time series x and y.

    Args:
        x (array-like)
        y (array-like)

    Returns:
        float: cov 
    """
    checkTS(x)
    checkTS(y)
    checkEqualLengths(x,y)
    return np.cov(x,y,ddof=1)[0,1]


def weighted_mean(data, err) -> float:
    """ 
    Compute the weighted mean. Here the weights are Gaussian error bars.
    See e.g. https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html.

    Args:
        data (array-like)
        errcov (array-like): error if 1d or cov if 2d. 

    Returns:
        float: weighted mean 
    """
    checkType(np.ndarray,data=data) 
    checkType(np.ndarray,err=err) 
    weights  = 1/err**2
    return np.sum( weights @ data )/np.sum(weights)


def weighted_variance(err) -> float:
    """ 
    Get variance of above weighted mean, when the weights are statistical errors. 

    Args:
        err (array-like)

    Returns:
        float: weighted variance 
    """
    checkType(np.ndarray,err=err) 
    weights = 1/err**2
    return 1/np.sum(weights)


def biased_sample_variance(data, err) -> float:
    """ 
    Compute the biased weighted sample variance, i.e. the biased variance of an 
    individual measurement and not the variance of the mean.

    Args:
        data (array-like)
        err (array-like)

    Returns:
        float: sample variance 
    """
    checkType(np.ndarray,data=data) 
    checkType(np.ndarray,err=err) 
    mean = weighted_mean(data, err)
    weights = 1/err**2
    V1 = np.sum(weights)
    return weights.dot((data - mean)**2) / V1


def unbiased_sample_variance(data, err) -> float:
    """ 
    Compute the unbiased weighted sample variance, i.e. the unbiased variance of an individual measurement and not
    the variance of the mean. Do not use this function if your weights are frequency weights. 
    """
    checkType(np.ndarray,data=data) 
    checkType(np.ndarray,err=err) 
    weights = 1/err**2
    V1 = np.sum(weights)
    V2 = np.sum(weights**2)
    return biased_sample_variance(data, weights) / ( 1 - (V2 / V1**2))


def unbiased_mean_variance(data, err) -> float:
    """ 
    Compute the unbiased variance of a weighted mean. Do not use this function if your weights are frequency
    weights. This is more like a systematic error. The absolute size of the weights does not matter. The error is
    constructed using the deviations of the individual data points. 
    """
    checkType(np.ndarray,data=data) 
    checkType(np.ndarray,err=err) 
    weights = 1/err**2
    V1 = np.sum(weights)
    V2 = np.sum(weights**2)
    return biased_sample_variance(data, weights) * V2 / ( V1**2 - V2)


def cov_to_cor(cov) -> np.ndarray:
    """ 
    Normalize covariance matrix to create correlation matrix.

    Args:
        cov (np.ndarray)

    Returns:
        np.ndarray: correlation matrix 
    """
    checkType(np.ndarray,cov=cov) 
    diagonal_sqrt = np.sqrt(np.diag(cov))
    return cov / np.outer(diagonal_sqrt, diagonal_sqrt)


def confidence_ellipse(x,y,ax,color='r',CI=None):
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
    checkType(np.ndarray,x=x)
    checkType(np.ndarray,y=y)
    checkEqualLengths(x,y)
    if CI is None:
        s = 1
    else:
        if CI<0 or CI>1:
            logger.TBRaise('Confidence limits are between 0 and 1.') 
        s = np.sqrt(-2*np.log(1-CI))
    data = np.vstack((x, y))
    cov  = np.cov(data)
    eigvals, eigvecs = np.linalg.eig(cov)
    maj_eigvec = eigvecs[:,np.argmax(eigvals)]
    theta = np.rad2deg( np.arctan2(maj_eigvec[1],maj_eigvec[0]) )
    a = s*np.sqrt(np.max(eigvals))
    b = s*np.sqrt(np.min(eigvals))
    ellipse = matplotlib.patches.Ellipse((std_mean(x), std_mean(y)), width = 2*a, height = 2*b, angle = theta,
                                          edgecolor=color, fc='None', zorder=FOREGROUND)
    ax.add_patch(ellipse)
    return a, b, theta 


def plot_correlation(mat,ax=plt):
    """ Plot correlation matrix as a heatmap.

    Args:
        mat (np.ndarray): correlation matrix
        ax (matplotlib ax object): Defaults to plt.
    """
    checkType(np.ndarray,mat=mat)
    heatmapColors = ['blue', 'white', 'red']
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom_cmap', heatmapColors)
    ax.imshow(mat, cmap=cmap, vmin=-1,vmax=1)
    ax.colorbar()


def forcePositiveSemidefinite(mat):
    """ Doctors a noisy correlation matrix mat to be positive semidefinite if it isn't already. 
    Uses algorithm of Rebonato and Jaeckel, DOI: 10.2139/ssrn.1969689
    
    Args:
        mat (np.ndarray)

    Returns:
        np.ndarray: positive semidefinite matrix 
    """
    eigvals, S = np.linalg.eig(mat)
    eigvals[eigvals<0] = 0
    sqrtLambda = np.diag(np.sqrt(eigvals))

    T = np.sum(S**2*eigvals,axis=1)
    sqrtT = np.diag(np.sqrt(1/T))
    B = sqrtT @ S @ sqrtLambda

    return B @ B.T


def dev_by_dist(data, axis=0, return_both_q=False, percentile=68):
    """ 
    Calculate the distance between the median and 68% quantiles. Returns the larger of the two 
    distances. This method is used sometimes to estimate error, for example in the bootstrap. 
    """
    data = np.asarray(data)
    median = np.nanmedian(data, axis)
    numb_data = data.shape[axis]
    idx_dn = max(int(np.floor((numb_data-1) / 2 - percentile/100 * (numb_data-1) / 2)), 0)
    idx_up = min(int(np.ceil((numb_data-1) / 2 + percentile/100 * (numb_data-1) / 2)), numb_data-1)
    #TODO: Need to handle the situation idx_dn=idx_up a bit carefully. the SRI project presents
    # some data with this issue, and you can troubleshoot that.
    sorted_data = np.sort(data - np.expand_dims(median, axis), axis=axis)
    q_l = np.take(sorted_data, idx_dn, axis)
    q_r = np.take(sorted_data, idx_up, axis)
    if return_both_q:
        return np.abs(q_l), np.abs(q_r)
    else:
        return np.max(np.stack((np.abs(q_l), np.abs(q_r)), axis=0), axis=0)


def error_prop(func, means, errors, grad=None, args=()):
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
    checkType(np.ndarray,means=means)
    checkType(np.ndarray,errors=errors)
    if isHigherDimensional(means):
        logger.TBRaise('means must be 1-d array.',frame=3)
    mean = func(means, *args)

    # Test if we got a covariance matrix
    if not isHigherDimensional(errors):
        errors = np.diag(errors**2)

    if type(mean) is tuple:
        logger.TBRaise("Tuples are not supported for error propagation")

    if grad is not None:
        grad = grad(means, *args)
    else:
        grad = diff_jac(means, func, args).T

    error = 0
    try:
        for i in range(len(grad)):
            for j in range(len(grad)):
                error += grad[i] * grad[j] * errors[i, j]
        error = np.sqrt(error)
    except TypeError:
        error += abs(grad * errors[0])
    return unvector(mean), unvector(error)


def error_prop_func(x, func, params, params_err, grad=None, args=()):
    """ 
    Propagate error in f(x;params,params_err). This needs its own special treatment, since
    the error propagation method on its own only propagates params_err to f(params,params_err).

    Args:
        x (array-like)
        func (func)
        params (array-like): Model parameters. 
        params_err (array-like): Error in model parameters. 
        grad (func, optional): Gradient function. Defaults to None.
        args (tuple, optional): Arguments of func. Defaults to ().
    """
    def wrap_func(m,a=()):
        if isinstance(a,tuple): 
            return func(x,m,*a)
        else:
            return func(x,m,a)
    if grad is None:
        wrap_grad = None
    else:
        def wrap_grad(m,a=()):
            if isinstance(a,tuple): 
                return grad(x,m,*a)
            else:
                return grad(x,m,a)
    return error_prop(wrap_func, params, params_err, wrap_grad, args)[1]


def gaudif(x1,e1,x2,e2) -> float:
    """ 
    Likelihood that difference between outcomes x1 and x2 is due to chance, assuming x1 and x2 are
    both drawn from a normal distribution with the same mean. A rule of thumb is that this is more
    appropriate when one estimated x1 and x2 using ~30 or more measurements.

    Args:
        x1 (float): mean 1 
        e1 (float): error 1
        x2 (float): mean 2
        e2 (float): error 2

    Returns:
        float: p-value 
    """
    if e1<0 or e2<0:
        logger.TBRaise('Error bars should be non-negative. Got',e1,e2)
    sigma = np.sqrt(e1**2 + e2**2)
    z     = abs(x1-x2)/(sigma * np.sqrt(2.))
    return 1.0 - sp.special.erf(z)


def studif(x1,e1,ndat1,x2,e2,ndat2) -> float:
    """ 
    Likelihood that difference between outcomes x1 and x2 is due to chance, assuming x1 and x2 are
    both drawn from a normal distribution with the same mean. A rule of thumb is that this is more
    appropriate when one estimated x1 and x2 using ~30 or fewer measurements.

    Args:
        x1  (float): mean 1 
        e1  (float): error 1
        ndat1 (int): number of measurements used to compute x1
        x2  (float): mean 2
        e2  (float): error 2
        ndat2 (int): number of measurements used to compute x2

    Returns:
        float: p-value 
    """
    if e1<0 or e2<0:
        logger.TBRaise('Error bars should be non-negative. Got',e1,e2)
    if ndat1<1 or ndat2 <1:
        logger.TBRaise('Need at least 2 data. Got',ndat1,ndat2)
    dof   = ndat1 + ndat2 - 2
    var12 = ndat1*e1**2
    var22 = ndat2*e2**2
    sigma = np.sqrt( (1/ndat1+1/ndat2)*( (ndat1-1)*var12 + (ndat2-1)*var22 )/dof )
    t     = abs(x1-x2)/sigma
    x     = dof/(dof+t**2)
    if t==0:
        return 1
    elif t>0:
        return sp.special.betainc(dof/2,1/2,x)
    else:
        return 2-sp.special.betainc(dof/2,1/2,x)


def goodnessOfFit(dof, chi2) -> float:
    """ 
    The q-value or goodness of fit.

    Args:
        dof (int): number of degrees of freedom 
        chi2 (float): the chi2 

    Returns:
        float: Q 
    """
    checkType(int,dof=dof)
    return sp.special.gammaincc(dof/2,chi2/2)


def plot_func(func, domain, params=(), args=(), func_err=None, params_err=(), 
              grad = None, swapXY=False, npoints=1000, **kwargs):
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
    checkDomain(domain)
    checkType(int,npoints=npoints)
    fill_param_dict(kwargs)
    kwargs['marker'] = None
    xmin = domain[0] 
    xmax = domain[1] 

    xdata = np.arange(xmin, xmax, (xmax - xmin) / npoints)
    try:
        ydata = func(xdata, params, *args)
    except TypeError:
        ydata = func(xdata, *params, *args)

    # Received an explicit error function:
    if func_err is not None:
        ydata_err = func_err(xdata, params, params_err, *args)
        if swapXY:
            return plot_fill(xdata, ydata, yedata=None, xedata=ydata_err, center=False, **kwargs)
        else:
            return plot_fill(xdata, ydata, ydata_err, center=False, **kwargs)

    # No explicit error function, but received a covariance matrix: 
    elif len(params_err) > 0:
        ydata_err = error_prop_func(xdata, func, params=params, params_err=params_err, grad=grad, args=args)
        if swapXY:
            return plot_fill(xdata, ydata, yedata=None, xedata=ydata_err, center=False, **kwargs)
        else:
            return plot_fill(xdata, ydata, ydata_err, center=False, **kwargs)

    # No errors at all: 
    else:
        if swapXY:
            return plot_lines(ydata, xdata, yedata=None, xedata=None, **kwargs)
        else:
            return plot_lines(xdata, ydata, yedata=None, xedata=None, **kwargs)


def save_func(func, domain, params=(), args=(), func_err=None, params_err=(), 
              grad = None, npoints=1000, header=None, filename='func.d', **kwargs):
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
    checkDomain(domain)
    checkType(int,npoints=npoints)
    checkType(str,filename=filename)
    xmin = domain[0] 
    xmax = domain[1] 

    createFilePath(filename)

    xdata = np.arange(xmin, xmax, (xmax - xmin) / npoints)
    try:
        ydata = func(xdata, params, *args)
    except TypeError:
        ydata = func(xdata, *params, *args)

    if func_err is not None:
        ydata_err = func_err(xdata, params, params_err, *args)
        writeTable(filename,xdata,ydata,ydata_err,header=header,**kwargs)

    elif len(params_err) > 0:
        ydata_err = error_prop_func(xdata, func, params=params, params_err=params_err, grad=grad, args=args)
        writeTable(filename,xdata,ydata,ydata_err,header=header,**kwargs)

    else:
        writeTable(filename,xdata,ydata,header=header,**kwargs)


def getModelWeights(IC) -> np.ndarray:
    """ 
    Convert information criteria IC to normalized probability weights.

    Args:
        IC (array-like): Array of information criteria 

    Returns:
        np.array: Probability weights 
    """
    checkType(np.ndarray,IC=IC)
    return normalize(np.exp(-0.5*IC))


def modelAverage(data,err,IC,return_syst=False):
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
    checkType(np.ndarray,data=data)
    checkType(np.ndarray,err=err)
    checkType(np.ndarray,IC=IC)
    checkEqualLengths(data,err,IC)
    pr = getModelWeights(IC)
    mean = np.sum(pr*data)
    var_stat = np.sum(pr*err**2)
    var_syst = np.sum(pr*data**2) - mean**2
    var = var_stat + var_syst 
    if return_syst:
        return mean, np.sqrt(var), np.sqrt(var_syst)
    else:
        return mean, np.sqrt(var)


def empiricalCDF(data):
    """ 
    Create the x and y coordinates needed to plot the empirical CDF
    of a 1-d set of data.

    Args:
        data (array-like): measurements 

    Returns:
        func: CDF 
    """
    checkTS(data)
    return sp.stats.ecdf(data).cdf.evaluate


def KSTest_2side(data1,data2) -> float:
    """ 
    2-sided Kolmogorov test. Gives back the likelihood that the observed difference between
    data1 and data2 are at least as extreme as suggested by the Kolmogorov statistic.

    Args:
        data1 (array-like)
        data2 (array-like)

    Returns:
        float: 1-p 
    """
    checkType(np.ndarray,data1=data1)
    checkType(np.ndarray,data2=data2)
    data1,data2 = toNumpy(data1,data2)
    return 1 - sp.stats.kstest(data1, data2).pvalue 


def KSTest_1side(data,cdf) -> float:
    """ 
    1-sided Kolmogorov test. Gives back the likelihood that the observed difference between
    data and cdf are at least as extreme as suggested by the Kolmogorov statistic.

    Args:
        data (array-like)
        cdf (function)

    Returns:
        float: 1-p 
    """
    checkType(np.ndarray,data=data)
    return 1 - sp.stats.kstest(data, cdf).pvalue


def binSeries(data,nbins) -> np.ndarray:
    """ 
    Take a time series and bin it. Bin 0 is the average over the first binsize elements,
    bin 1 the average over the next binsize elements, and so on.

    Args:
        data (array-like)
        nbins (int)

    Returns:
        np.ndarray: Binned data
    """
    checkTS(data)
    checkType(int,nbins=nbins)
    ndat=len(data)
    binsize=int(ndat/nbins)
    reshaped_data = data[:binsize * nbins].reshape(nbins, binsize)
    return np.apply_along_axis(std_mean, 1, reshaped_data)


def symmetrizeError(lo,hi,central,method='conservative') -> float:
    checkType(str,method=method)
    checkType("real",lo=lo)
    checkType("real",hi=hi)
    checkType("real",central=central)
    if (lo<=0) or (hi<=0):
        logger.TBRaise('Uncertainties must be positive. Got lo, hi =',lo,hi)
    if method=='conservative':
        err  = np.max([lo,hi])
        mean = central
    # FLAG method from 10.1140/epjc/s10052-022-10536-1
    elif method=='FLAG':
        if lo>hi:
            mean = central + (lo-hi)/4
            err  = (hi + 3*lo)/4
        else:
            mean = central + (hi-lo)/4
            err  = (lo + 3*hi)/4
    else:
        logger.TBRaise('Unrecognized method',method)
    return mean, err
