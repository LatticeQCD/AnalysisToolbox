#
# statistics.py
#
# H. Sandmeyer, D. Clarke
#
# A collection of basic methods for statistical analysis. The following methods have been adapted from software of
# Bernd Berg, Markov Chain Monte Carlo and their Statistical Analysis, World Scientific, 2004, ISBN=978-981-3106-37-6:
#     gaudif, studif
#


import numpy as np
from scipy.linalg import inv
from scipy.special import betainc, erf
from latqcdtools.math.num_deriv import diff_jac 
from latqcdtools.math.math import logDet
from latqcdtools.base.plotting import fill_param_dict, plot_fill, plot_lines
from latqcdtools.base.utilities import isHigherDimensional
from latqcdtools.base.cleanData import clipRange
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger


def reduce_tuple(func):
    def func_wrapper(data, *args, **kwargs):
        if type(data[0]) is tuple:
            retvalue = ()
            for i in range(len(data[0])):
                obj_array = []
                for k in range(len(data)):
                    obj_array.append(data[k][i])
                retvalue += (func(obj_array, *args, **kwargs),)
            return retvalue
        else:
            return func(data, *args, **kwargs)
    return func_wrapper


def meanArgWrapper(func,used_data,args):
    if isinstance(args, dict):
        return func(used_data, **args)
    else:
        return func(used_data, *args)


@reduce_tuple
def std_median(data, axis = 0):
    """ Compute the median. The default behavior of numpy is to flatten the data, flagged by axis=None. This can be
    inconvenient, for example in the bootstrap and jackknife routines. It is also inconvenient if you have e.g. an
    array of matrices, and you want the median matrix rather than the median element. """
    return np.median(data, axis)


@reduce_tuple
def std_mean(data, axis = 0):
    """ Compute the mean. """
    return np.mean(data, axis)


@reduce_tuple
def std_var(data, axis = 0):
    """ Calculate unbiased (ddof = 1) estimator for the variance. """
    data = np.asarray(data)
    return np.var(data, axis = axis, ddof = 1)


@reduce_tuple
def std_dev(data, axis = 0):
    """ Calculate unbiased (ddof = 1) estimator for the standard deviation. """
    data = np.asarray(data)
    return np.std(data, axis = axis, ddof = 1)


@reduce_tuple
def std_err(data, axis = 0):
    """ Standard deviation of the sample mean, according the the CLT. """
    data = np.asarray(data)
    return std_dev(data, axis) / np.sqrt(data.shape[axis])


def expandArgs(func, x, params=(), args=()):
    """ In general we distinguish between parameters and arguments. Parameters should be passed
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
    """ Some methods require that you do something over an interval, which we refer to in this module as
    a 'domain'. This checks the domain makes sense.

    Args:
        domain (tuple)
    """
    checkType(domain,tuple)
    if len(domain) != 2:
        logger.TBError('A domain is a tuple of the form (xmin,xmax) specifying the interval [xmin,xmax].')
    if domain[0]>=domain[1]:
        logger.TBError('Must have domain[1]>domain[0].')


def checkPrior(prior,prior_err):
    """ Make sure prior and prior_err status are compatible. """
    if prior is None and prior_err is not None:
        logger.TBError('prior = None, prior_err != None')
    if prior_err is None and prior is not None:
        logger.TBError('prior != None, prior_err = None')


def checkTS(ts):
    """ Some methods require 1-d time series. This checks that the type, dimensionality,
    and length are appropriate.

    Args:
        ts (array-like): time series 
    """
    checkType(ts,'array')
    if isHigherDimensional(ts):
        logger.TBError('Expected 1-d time series.')
    if len(ts) < 2:
        logger.TBError('Time series needs at least two measurements.')


def countParams(func,params):
    """ If we have a function, return length of params. Else, we must have a spline. """
    nparams = len(params)
    if nparams == 0:
        try:
            # LSQUnivariateSpline
            nparams = len(func.get_coeffs())
            # CubicSpline
        except AttributeError:
            nparams = len(func.c)
    return nparams


def DOF(ndat,nparam,priorsigma):
    """ Compute the number of degrees of freedom. Depends on whether you use priors. Any input priors are taken as
    initial guesses for the fit algorithm. If you would like parameter in the prior array to be treated as a 
    starting guess only, and not as a Bayesian prior, set its corresponding error to np.inf. Hence when there
    are priors, the number of degrees of freedom equals the number of ydata, less the number of finite prior errors. """
    if priorsigma is not None:
        dof = ndat - np.count_nonzero(priorsigma == np.inf)
    else:
        dof = ndat - nparam 
    if dof < 0:
        logger.TBError('Fewer data than fit parameters!')
    logger.debug('Computed d.o.f. =',dof)
    return dof


def chisquare(xdata,ydata,cov,func,args=(),params=(),prior=None,prior_err=None):
    """ Calculate chi^2.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): fit parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        prior_err (array-like, optional): Bayesian prior errors. Defaults to None.

    Returns:
        float: chi^2 
    """
    checkPrior(prior,prior_err)
    y    = expandArgs(func,xdata,params,args)
    cor  = norm_cov(cov)
    diff = ( ydata - y )/np.sqrt( np.diag(cov) )
    res  = diff.dot( inv(cor).dot(diff) )
    if prior is not None:
        res += np.sum((np.array(params) - prior)**2 / prior_err**2)
    return res


def logGBF(xdata, ydata, cov, func, args=(), params=(), prior=None, prior_err=None):
    """ log P(data|model). This quantity is useful for comparing fits of the same data to different models that
    have different priors and/or fit functions. The model with the largest logGBF is the one preferred by the data.
    Differences in logGBF smaller than 1 are not very significant. Gaussian statistics are assumed.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): fit parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        prior_err (array-like, optional): Bayesian prior errors. Defaults to None.

    Returns:
        float: log( Gaussian Bayes factor ) 
    """
    chi2    = chisquare(xdata, ydata, cov, func, args, params, prior, prior_err)
    nparams = countParams(func,params)
    dof     = DOF(len(ydata),nparams,prior_err)
    logger.debug('chi^2 =',chi2)
    logger.debug('nparams =',nparams)
    logger.debug('dof =',dof)
    if prior is None:
        return 0.5*( - logDet(cov) - chi2 - dof*np.log(2*np.pi) )
    else:
        return 0.5*( - logDet(cov) - chi2 - dof*np.log(2*np.pi) + logDet(np.diag(clipRange(prior_err)**2)) )


def AIC(xdata, ydata, cov, func, args=(), params=(), prior=None, prior_err=None):
    """ The Akaike information criterion (AIC) is a measure of how well a fit performs. It builds on the likelihood
    function by including a penalty for each d.o.f. This is useful in a context where you have multiple models to
    choose from,and hence different numbers of d.o.f. possible. It's also useful when you are worried about
    overfitting. The preferred model minimizes the AIC.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): fit parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        prior_err (array-like, optional): Bayesian prior errors. Defaults to None.

    Returns:
        float: AIC
    """
    nparams    = countParams(func,params)
    likelihood = logGBF(xdata, ydata, cov, func, args, params, prior, prior_err)
    return 2*nparams - 2*likelihood


def AICc(xdata, ydata, cov, func, args=(), params=(), prior=None, prior_err=None):
    """ Corrected AIC (AICc). When the sample size is smaller, it increases the chance AIC will select a model with too
    many parameters. The AICc tries to further correct for this. In the limit that the number of data points goes to
    infinity, one recovers the AIC.

    Args:
        xdata (array-like)
        ydata (array-like)
        cov (array-like): covariance matrix 
        func (func)
        args (tuple, optional): arguments to func. Defaults to ().
        params (tuple, optional): fit parameters. Defaults to ().
        prior (array-like, optional): Bayesian priors. Defaults to None.
        prior_err (array-like, optional): Bayesian prior errors. Defaults to None.

    Returns:
        float: corrected AIC 
    """
    nparams = countParams(func,params)
    ndat    = len(ydata)
    aic     = AIC(xdata, ydata, cov, func, args, params, prior, prior_err)
    return aic + 2*(nparams**2+nparams)/(ndat-nparams+1)


def pearson(x,y):
    """ Get the Pearson correlation coefficient between the time series x and y.

    Args:
        x (array-like): _description_
        y (array-like): _description_

    Returns:
        float: R 
    """
    checkTS(x)
    checkTS(y)
    return np.corrcoef(x,y)[0,1]


def weighted_mean(data, err):
    """ Compute the weighted mean. Here the weights are Gaussian error bars.
    See e.g. https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html.

    Args:
        data (array-like)
        errcov (array-like): error if 1d or cov if 2d. 

    Returns:
        float: weighted mean 
    """
    data = np.array(data)
    weights  = 1/np.array(err)**2
    return np.sum( weights @ data )/np.sum(weights)


def weighted_variance(err):
    """ Get variance of above weighted mean, when the weights are statistical errors. 

    Parameters
    ----------

    err: array_like
        The errors of the data points.
    """
    weights = 1/np.array(err)**2
    return 1/np.sum(weights)


def biased_sample_variance(data, err):
    """ Compute the biased weighted sample variance, i.e. the biased variance of an individual measurement and not the
    variance of the mean. """
    mean = weighted_mean(data, err)
    weights = 1/np.array(err)**2
    V1 = np.sum(weights)
    return weights.dot((data - mean)**2) / V1


def unbiased_sample_variance(data, err):
    """ Compute the unbiased weighted sample variance, i.e. the unbiased variance of an individual measurement and not
    the variance of the mean. Do not use this function if your weights are frequency weights. """
    weights = 1/np.array(err)**2
    V1 = np.sum(weights)
    V2 = np.sum(weights**2)
    return biased_sample_variance(data, weights) / ( 1 - (V2 / V1**2))


def unbiased_mean_variance(data, err):
    """ Compute the unbiased variance of a weighted mean. Do not use this function if your weights are frequency
    weights. This is more like a systematic error. The absolute size of the weights does not matter. The error is
    constructed using the deviations of the individual data points. """
    data    = np.asarray(data)
    weights = 1/np.array(err)**2
    V1 = np.sum(weights)
    V2 = np.sum(weights**2)
    return biased_sample_variance(data, weights) * V2 / ( V1**2 - V2)


def norm_cov(cov):
    """ Normalize a covariance matrix to create the correlation matrix. """
    res = np.zeros((len(cov), len(cov[0])))
    for i in range(len(cov)):
        for j in range(len(cov[0])):
            res[i][j] = cov[i][j] / np.sqrt( cov[j][j] * cov[i][i] )
    return np.array(res)


def cut_eig(corr, threshold):
    """ Cut eigenvalues of the correlation matrix. If they are smaller than the threshold, replace them with the
    threshold. When needed, this replaces a small eigenvalue by a larger, small eigenvalue, which has the effect of
    slightly overestimating the errors. The alternative would be to ignore them, in which case the program would
    crash because the matrix is singular, or to discard them, which is like setting the variance to infinity.
    This procedure is more accurate than the latter option. """
    vals, vecs = np.linalg.eig(corr)
    for i, value in enumerate(vals):
        if value < threshold:
            logger.details('Set small eigenvalue',value,'from correlation matrix to threshold',threshold)
            vals[i] = threshold
    return vecs.dot( np.diag(vals).dot( vecs.transpose() ) )


@reduce_tuple
def dev_by_dist(data, axis=0, return_both_q=False, percentile=68):
    """ Calculate the distance between the median and 68% quantiles. Returns the larger of the two distances. This
    method is used sometimes to estimate error, for example in the bootstrap. """
    data = np.asarray(data)
    median = np.nanmedian(data, axis)
    numb_data = data.shape[axis]
    idx_dn = max(int(np.floor((numb_data-1) / 2 - percentile/100 * (numb_data-1) / 2)), 0)
    idx_up = min(int(np.ceil((numb_data-1) / 2 + percentile/100 * (numb_data-1) / 2)), numb_data-1)
    sorted_data = np.sort(data - np.expand_dims(median, axis), axis=axis)
    q_l = np.take(sorted_data, idx_dn, axis)
    q_r = np.take(sorted_data, idx_up, axis)
    if return_both_q:
        return np.abs(q_l), np.abs(q_r)
    else:
        return np.max(np.stack((np.abs(q_l), np.abs(q_r)), axis=0), axis=0)


def error_prop(func, means, errors, grad=None, args=()):
    """ Use error propagation to propagate some errors through function func. The function should have the form
        func( data ), where data is your array of input variables. """
    errors = np.asarray(errors)
    means  = np.asarray(means)
    mean   = func(means, *args)

    # Test if we got a covariance matrix
    if not isHigherDimensional(errors):
        errors = np.diag(errors**2)

    if type(mean) is tuple:
        raise TypeError("Tuples are not supported for error propagation")

    if grad is not None:
        grad = grad(means, *args)
    else:
        grad = diff_jac(means, func, args).transpose()
    error = 0
    try:
        for i in range(len(grad)):
            for j in range(len(grad)):
                error += grad[i] * grad[j] * errors[i, j]
        error = np.sqrt(error)
    except TypeError:
        error += abs(grad * errors[0])
    return mean, error


def error_prop_func(x, func, params, params_err, grad=None, args=()):
    """ Propagate error in f(x;params,params_err). This needs its own special treatment, since
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


#def error_budget(x, func, means, errors, grad=None, args=()):
#    if isHigherDimensional(errors):
#        logger.TBError('Covariance matrix not yet supported.')
#    # Strategy is to see separately contribution from each error
#    sigma_f0 = error_prop_func(x,func,means,errors,grad=grad,args=args)
#    sigma_fi = []
#    for i in range(len(errors)):
#        err_contrib    = np.zeros(len(errors))
#        err_contrib[i] = errors[i]
#        sigma_fi.append(error_prop_func(x,func,means,err_contrib,grad=grad,args=args))
#    sigma_fi = np.array(sigma_fi)
#    perc_fi = sigma_fi/sigma_f0
#    print(np.sum(perc_fi))
#    for i in range(len(errors)):
#        logger.info('x_'+str(i),':',round(perc_fi[i],4))


def gaudif(x1,e1,x2,e2):
    """ Likelihood that difference between outcomes x1 and x2 is due to chance, assuming x1 and x2 are
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
        logger.TBError('Error bars should be non-negative. Got',e1,e2)
    sigma = np.sqrt(e1**2 + e2**2)
    z     = abs(x1-x2)/(sigma * np.sqrt(2.))
    return 1.0 - erf(z)


def studif(x1,e1,ndat1,x2,e2,ndat2):
    """Likelihood that difference between outcomes x1 and x2 is due to chance, assuming x1 and x2 are
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
        logger.TBError('Error bars should be non-negative. Got',e1,e2)
    if ndat1<1 or ndat2 <1:
        logger.TBError('Need at least 2 data. Got',ndat1,ndat2)
    dof   = ndat1 + ndat2 - 2
    var12 = ndat1*e1**2
    var22 = ndat2*e2**2
    sigma = np.sqrt( (1/ndat1+1/ndat2)*( (ndat1-1)*var12 + (ndat2-1)*var22 )/dof )
    t     = abs(x1-x2)/sigma
    x     = dof/(dof+t**2)
    if t==0:
        return 1
    elif t>0:
        return betainc(dof/2,1/2,x)
    else:
        return 2-betainc(dof/2,1/2,x)


def plot_func(func, domain, params=(), args=(), func_err=None, params_err=(), 
              grad = None, swapXY=False, npoints=1000, **kwargs ):
    """ Plot a function along with its error bands.

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
    fill_param_dict(kwargs)
    kwargs['marker'] = None
    checkDomain(domain)
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
            return plot_fill(xdata, ydata, yedata=None, xedata=ydata_err, **kwargs)
        else:
            return plot_fill(xdata, ydata, ydata_err, **kwargs)

    # No explicit error function, but received a covariance matrix: 
    elif len(params_err) > 0:
        ydata_err = error_prop_func(xdata, func, params=params, params_err=params_err, grad=grad, args=args)
        if swapXY:
            return plot_fill(xdata, ydata, yedata=None, xedata=ydata_err, **kwargs)
        else:
            return plot_fill(xdata, ydata, ydata_err, **kwargs)

    # No errors at all: 
    else:
        if swapXY:
            return plot_lines(ydata, xdata, yedata=None, xedata=None, **kwargs)
        else:
            return plot_lines(xdata, ydata, yedata=None, xedata=None, **kwargs)
