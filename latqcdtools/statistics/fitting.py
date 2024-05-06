#
# fitting.py
#
# H. Sandmeyer, D. Clarke
#
# Fitter class, along with functions for fitting data. The goal is to be able to carry out simple fits, trying many
# fitting algorithms. Bayesian priors are supported. You can judge the quality of your fit by quantities like
# chi^2/d.o.f. and logGBF. 
#

import numpy as np
from scipy.optimize import curve_fit
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkEqualLengths
from latqcdtools.base.speedify import DEFAULTTHREADS, parallel_function_eval
from latqcdtools.base.plotting import plot_dots, plot_bar, plt
from latqcdtools.base.readWrite import writeTable
from latqcdtools.base.utilities import envector, isHigherDimensional, toNumpy, createFilePath
from latqcdtools.math.math import invert, regulate, checkSquare, isSymmetric, isPositiveSemidefinite
from latqcdtools.math.optimize import minimize
from latqcdtools.math.num_deriv import diff_jac, diff_fit_grad
from latqcdtools.statistics.statistics import plot_func, error_prop_func, cov_to_cor, chisquare, logGBF, DOF, \
    expandArgs, checkDomain, BAIC, AIC, AICc, checkPrior, goodnessOfFit


# Allowed keys for the constructor
_allowed_keys = ['grad', 'args', 'grad_args', 'tol', 'use_diff', 'error_strat',
                'norm_err_chi2', 'derive_chisq', 'svdcut', 'test_tol', 'max_fev', 'nproc']

# All possible algorithms.
all_algs = ["curve_fit", "L-BFGS-B", "TNC", "Powell", "Nelder-Mead", "COBYLA", "SLSQP", "CG","dogleg", "trust-ncg"]

# Standard algorithms for the minimization. 
std_algs = ["curve_fit", "TNC", "Powell", "Nelder-Mead"]

# Fast algorithms that work with priors. 
bayes_algs = ["TNC", "Powell", "Nelder-Mead"]


class Fitter:
    """ The :class:`Fitter`, contains all information necessary for fitting: The data, the function to be fitted, and
    optional the data for the errors. There are different minimization algorithms available. Many of them need the
    gradient or hessian of the chisquare. One way is to set the derivatives of the fitting function from outside.
    The derivatives of the actual chisquare are then computed via error propagation. Another way is to use numerical
    derivatives.

    There are two ways to compute the derivatives of the chisqare numerically. Either compute the
    numerical derivative of the whole chisquare (error_strat='hessian') or compute the derivatives of the fitting
    function and use error propagation (error_strat='propagation'). The latter is the default case."""

    def __init__(self, func, xdata, ydata, edata = None, **kwargs):
        """
        Parameters
        ----------
        func : callable
            Function to be fitted. Must have form 
                func(x, params, *args)
        xdata : array_like
            xdata used for fit. These data may be higher dimensional. This may be the case when our fit functions needs
            more than one parameter. However, the number of elements in the first axis has to be equal to the number of
            elements in ydata.
        ydata : array_like
            ydata used for fit.
        edata : array_like, optional, default: None
            Data for the error. Either pass an 1D array of errors of a full covariance matrix. In case of errors, the
            errors are interpreted as edata = sqrt(variance). For the case of the covariance matrix no root has to be
            taken: variance = diag(edata).
        grad : callable, optional, default: None
            gradient of the fit function.
        args : array_like, optional, default: ()
            Optional arguments that shall be passed to func and that should not be fitted.
        grad_args : array_like, optional, default: None
            Optional parameter for the gradient. If set to None the arguments for the function are used (args).
        tol : float, optional, default: 1e-12
            Tolerance for the minimization.
        svdcut : float, optional, default: 1e-12
            Apply this cut to the eigenvalues of the correlation matrix.
        max_fev : int, optional, default: 10000
            Maximum number of iterations / function evaluations.
        use_diff : bool, optional, default: True
            In case of numerical derivative use the difference quotient for approximation.
        norm_err_chi2 : bool, optional, default: False 
            Multiply errors with chi**2/d.o.f. Some people like to do this. 
        derive_chisq : bool, optional, default: False
            In case of numerical derivative, apply the derivative to the whole chisquare instead of the function.
        nproc : int, optional, default: DEFAULTTHREADS
            If you want you can accelerate the fits using nprocs threads.
        """

        diff = set(set(kwargs.keys()) - set(_allowed_keys))
        if len(diff) != 0:
            logger.TBError("Illegal argument(s) to fitter", *diff)

        # Some attributes that are set in functions other than __init__.
        self._grad = None
        self.grad  = None
        self._pcov = None

        # Store data
        checkEqualLengths(xdata,ydata)
        self._xdata, self._ydata = toNumpy(xdata, ydata) 

        # These attributes are described in the above doccumentation. If they aren't specified in the keyword
        # arguments when the Fitter is initialized, they take the default value shown here. 
        self._use_diff      = kwargs.get('use_diff', True)
        self._derive_chisq  = kwargs.get('derive_chisq', False)
        self._tol           = kwargs.get('tol', 1e-10)
        self._test_tol      = kwargs.get('test_tol', 1e-10)
        self._svdcut        = kwargs.get('svdcut', 1e-12)
        self._max_fev       = kwargs.get('max_fev', None)
        self._norm_err_chi2 = kwargs.get('norm_err_chi2', False)
        self._args          = kwargs.get('args', ())
        self._grad_args     = kwargs.get('grad_args', None)
        self._errorAlg      = kwargs.get('error_strat', 'propagation')
        self._nproc         = kwargs.get('nproc', 1)
        logger.debug('Initialize fitter with:')
        logger.debug('  use_diff:',self._use_diff) 
        logger.debug('  derive_chisq:',self._derive_chisq)
        logger.debug('  tol:',self._tol)
        logger.debug('  test_tol:',self._test_tol)
        logger.debug('  max_fev:',self._max_fev)
        logger.debug('  norm_err_chi2:',self._norm_err_chi2)
        logger.debug('  args:',self._args)
        logger.debug('  grad_args:',self._grad_args)
        logger.debug('  errorAlg:',self._errorAlg)
        logger.debug('  nproc:',self._nproc) 

        if (self._errorAlg=='hessian') and self._derive_chisq is True:
            self._errorAlg='propagation'
            logger.warn('hessian strategy not yet compatible with derive_chisq; switching to propagation.')

        if self._grad_args is None:
            self._grad_args = self._args

        if type(self._max_fev) is int:
            tmp_fev = self._max_fev
            self._max_fev = dict()
            for alg in all_algs:
                self._max_fev[alg] = tmp_fev
        if self._max_fev is None:
            self._max_fev = {
                "curve_fit"  : 50000,
                "L-BFGS-B"   : 15000,
                "TNC"        : 5000 ,
                "Powell"     : 10000,
                "Nelder-Mead": 5000 ,
                "COBYLA"     : 15000,
                "SLSQP"      : 15000,
                "CG"         : 15000,
                "dogleg"     : 15000,
                "trust-ncg"  : 15000
                }

        # Initialize func. This is also done in set_func, but we need it before that
        self._func = func

        # The final parameters and covariance matrices will be saved in this. 
        self._saved_params = None 
        self._saved_pcov   = None 

        self.set_func(func, kwargs.get('grad', None)) 

        # For Bayesian fits. If these are used, they are set in general_fit, and then utilized elsewhere.
        self._priorval   = None 
        self._priorsigma = None 

        # Check if we have the covariance matrix available. 
        if edata is not None:
            edata = np.asarray(edata, dtype = float)
            if isHigherDimensional(edata): 
                self._cov = edata
            else: 
                self._cov = np.diag(np.array(edata)**2)
        else: 
            # Initialize everything to one if we don't get error information. 
            self._cov = np.diag(np.ones(len(self._ydata)))

        checkSquare(self._cov)
        if not isSymmetric(self._cov):
            logger.TBError('Covariance matrix is not symmetric.')
        if not isPositiveSemidefinite(self._cov):
            logger.TBError('Covariance matrix is not positive semidefinite.')

        # Compute a correlation matrix. Note that if the incoming data are highly correlated, then
        # the covariance matrix is likely to be ill-conditioned. This results in numerical stability
        # issues, most obviously when trying to invert the matrix, but also when trying to do things
        # like compute its jacobian. Therefore regulate it. 
        self._cov     = regulate( self._cov, svdcut=self._svdcut )
        self._cor     = cov_to_cor( self._cov )
        self._inv_cor = invert( self._cor, 'scipy' )
        self._inv_cov = invert( self._cov, 'scipy' )


    def __repr__(self) -> str:
        return "Fitter"


    def set_func(self, func, grad = None, args = None, grad_args = None, hess_args = None):
        """ Set fitting function, gradient, and their arguments. Also initialize self.func,
        and self.grad. These point to the actual wrappers which are used in the fit. In
        case of provided gradient, this will be wrap_grad. In case of numerical derivatives, 
        this will be num_grad.

        Parameters
        ----------
        func : callable
            Function to be fitted.
        grad : callable, optional, default: None
            gradient of the fit function.
        args : array_like, optional, default: ()

        grad_args : array_like, optional, default: None
            Optional parameter for the gradient. If set to None the arguments for the function
            are used (args).
        """

        # Direct storage of the user functions.
        self._func = func
        self._grad = grad

        # Later we only access self.func and self.grad. These are wrappers around the user function or the
        # numerical derivatives. The code below chooses the right wrappers. When they are set to None, None will
        # ultimately be passed to scipy.optimize.minimize, which indicates that minimize should estimate these
        # quantities with numerical derivatives.
        if self._grad is None:
            if self._use_diff and not self._derive_chisq:
                self.grad = self.num_grad
            else:
                self.grad = None
                self._derive_chisq = True
        else:
            self.grad = self.wrap_grad

        if args is not None:
            self._args = args

        if grad_args is not None:
            self._grad_args = grad_args
        elif self._args is not None:
            self._grad_args = self._args


    def num_grad(self, x, params):
        return diff_fit_grad(x, params, self._func, self._args)


    def wrap_grad(self, x, params):
        return expandArgs(self._grad,x,params,self._grad_args)
    def wrap_func(self, x, params):
        return expandArgs(self._func,x,params,self._args)


    def fit_ansatz_array(self, params):
        """ Return the array of the fit ansatz values at each position in self._xdata. """
        params = np.array(params)
        return self.wrap_func(self._xdata, params)


    def jacobian_fit_ansatz_array(self, params):
        """ If f is the fit function, return the array df / dp_i evaluated at each _fit_xdata. """
        params = np.array(params)
        # Note that the function is passed in secret here.
        return self.grad(self._xdata, params).T


    def calc_chisquare(self, params):
        """ Compute chi^2. This is the function that will be minimized. """
        return chisquare(self._xdata, self._ydata, self._cov, self._func, self._args, params, 
                         prior=self._priorval,priorsigma=self._priorsigma)


    def grad_chisquare(self, params):
        """ Compute the gradient of the chisquare, see e.g. eq. (A3) in 10.1103/PhysRevD.90.054506. 
        We assume priors are not correlated with data. """
        jac  = self.jacobian_fit_ansatz_array(params).T     # df/dp
        diff = self.fit_ansatz_array(params) - self._ydata  # Dy
        res  = 2 * jac @ self._inv_cov @ diff
        if self._priorsigma is not None:
            res += np.sum( 2*( np.array(params) - self._priorval ) / self._priorsigma**2 )
        return res


    def hess_chisquare(self, params):
        """ Compute the Hessian of the chisquare. Terms of O(f(p)-y) are discarded, which is reasonable when
        the fit is relatively good. This also ensures non-negative eigenvalues for pcov and saves computational
        time computing the Hessian. See e.g. eq. (A4) in 10.1103/PhysRevD.90.054506. We assume priors are not 
        correlated with data. """
        jac = self.jacobian_fit_ansatz_array(params).T
        res = 2*( jac @ self._inv_cov @ jac.T )
        if self._priorsigma is not None:
            res += np.sum( 2/self._priorsigma**2 )
        return res


    def _num_func_jacobian(self, params):
        """ For the error computation we need the Jacobian of the array of function values. If self._derive_chisq is True,
        we cannot use self.grad_fit_ansatz_array. In that case, the Jacobian is calculated using this function. """
        return diff_jac(params, self.fit_ansatz_array)


    def minimize_chi2(self, start_params, algorithm):
        """ Minimize the chi^2 using the scipy minimize routine is used.

        Parameters
        ----------
        start_params: array_like
            The start parameters that are used for the minimization.
        algorithm: string
            The algorithm that will be used.

        Returns
        -------
        params:
            Array of the parameters that minimize the chisquare.
        chi2:
            The minimum value of the chisquare.
        """

        if self.grad is not None:
            jac_func=self.grad_chisquare
            hess_func=self.hess_chisquare
        else:
            jac_func=None
            hess_func=None

        if algorithm == "curve_fit":

            if self._priorsigma is not None:
                logger.TBError('The curve_fit algorithm is not yet able to handle priors.')

            # This ensures that however func is used, whatever is passed to him as arguments will be captured as a
            # tuple, which is then plugged into wrap_func.
            def func(x, *p):
                return self.wrap_func(x, p)

            cov = self._cov
            # If no gradient has been provided by the user, it is probably better to use the numerical derivative from
            # curve_fit instead of our own.
            if self._grad is not None:
                def grad(x, *p):
                    return np.array(self.grad(x, p)).T
            else:
                grad = None

            params, _ = curve_fit(func, self._xdata, self._ydata, sigma = cov, p0 = start_params, jac = grad,
                                  ftol = self._tol, maxfev = self._max_fev["curve_fit"])

        else:
            params = minimize(self.calc_chisquare, jac_func, hess_func, start_params, self._tol,
                              self._max_fev[algorithm], algorithm = algorithm)

        if any(np.isnan(params)) or any(np.isinf(params)):
            logger.TBFail(algorithm,": Fit result is inf or nan!")
            raise ValueError

        chi2 = self.calc_chisquare(params)

        return params, chi2


    def compute_err(self, params, chi2, algorithm):
        """ Compute the covariance matrix of the fit parameters. If no errors have been provided, they are assumed to
        be one. We get the fit variances from the diagonal elements of the covariance matrix.

        Parameters
        ----------
        params: array_like
            Parameters for which the errors should be computed.
        chi2: float
            The chisquare for those parameters.
        algorithm: string
            The algorithm used to compute these errors (Necessary for strings in exceptions).

        Returns
        ------
        pcov:
            Correlation matrix of the fit parameters.
        fit_errors:
            Diagonal of pcov.
        """

        dof = DOF(len(self._ydata),len(params),self._priorsigma)    

        logger.debug('using error alg',self._errorAlg)
        if dof <= 0:
            pcov = np.full((len(params),len(params)), np.nan)
        else:
            if self._errorAlg=='propagation':
                pcov = self.pcov_error_prop(params,algorithm)
            elif self._errorAlg=='hessian':
                pcov = self.pcov_hessian(params,algorithm)
            else:
                logger.TBError('Unknown fitting algorithm', self._errorAlg)

        # Sometimes people like to rescale the parameter covariance matrix by the chi^2/dof. This tries to
        # take the fit quality in account into the error directly, and to my understanding it is what
        # gnuplot does. In a physics context, this procedure seems to be unjustified, and it doesn't come 
        # naturally out of the mathematics, so this is not done by default.
        if self._norm_err_chi2:
            pcov *= chi2/dof
        fit_errors = np.sqrt(np.diag(pcov))

        self._pcov = pcov

        return pcov, fit_errors


    def pcov_error_prop(self, params, algorithm):
        """ Compute the parameter's covariance matrix through error propagation, i.e. pcov = (J^t * C^-1 * J)^-1, where
        J is the Jacobian of the fit function and C is the covariance matrix of the data points. """
        if self.grad is None:
            jac = self._num_func_jacobian(params)
        else:
            jac = self.jacobian_fit_ansatz_array(params)

        jej = np.matrix( jac.T @ self._inv_cov @ jac )

        try:
            pcov = invert(jej,'scipy')

            # Test that the inversion went relatively OK
            test = pcov @ jej
            if abs(np.sum(test) - np.sum(np.diag(test))) > self._test_tol:
                logger.warn(algorithm,"pcov @ jej not close to id")

            if np.min(np.diag(pcov)) < 0:
                logger.TBFail(algorithm + ": Negative entries for the variance!")
                raise ValueError

        except Exception as e:
            logger.TBFail(algorithm,"hit exception",e)
            raise e 

        return pcov 


    def pcov_hessian(self, params, algorithm):
        """ Obtain the parameter's covariance matrix by inverting the hessian of the chi^2. """
        pcov = invert(self.hess_chisquare(params)/2.)
        if np.min(np.diag(pcov)) < 0:
            logger.TBFail(algorithm + ": Negative entries for the variance!")
            raise ValueError
        return pcov


    def _general_fit(self, algorithm="curve_fit"): 
        """ Perform fit. No new fit data are generated. """
        params, chi2 = self.minimize_chi2(self._saved_params, algorithm)
        pcov, fit_errors = self.compute_err(params, chi2, algorithm)
        return params, fit_errors, chi2, pcov


    def _tryAlgorithm(self,algorithm):
        """ Wrapper that collects general fit results. Allows for parallelization. """
        try:
            params, fit_errors, chi2, pcov = self._general_fit(algorithm) 
            logger.details(algorithm, "successful. Chi^2 = ", chi2)
            return params, fit_errors, chi2, pcov, None
        except Exception as e:
            logger.details(algorithm, "failed with:", e)
            return None, None, np.inf, None, e


    def _autoDomain(self,domain):
        """ Set domain automatically if domain=None. """
        if domain is None:
            domain=(np.min(self._xdata),np.max(self._xdata))
        checkDomain(domain)
        return domain


    def try_fit(self, algorithms = std_algs, start_params = None, priorval = None, priorsigma = None, detailedInfo = False):
        """ Perform the fit. This is what you should usually call. Try different algorithms and choose the one with the
        smallest chi^2. By default this method does a standard statistical fit. One can also include priors to obtain
        posteriors using Bayesian statistics. A well known summary of the latter strategy in the context of lattice QCD
        is given here: https://arxiv.org/abs/hep-lat/0110175. If there are priors, calculate logGBF and return that too.

        Parameters
        ----------
        algorithms: string or None, optional, default = None
            List of strings with the algorithms that can be used. Possible values are:
                "L-BFGS-B", "TNC", "Powell" ,"Nelder-Mead", "COBYLA", "SLSQP", "CG", "BFGS", "dogleg", "trust-ncg".
            The latter 4 usually don't work. When provided None, the first 7 algorithms are used.
        start_params: array_like, optional, default = None
            The start parameters for the fit.
        priorval: array_like, optional, default = None
            For constrained fits. Prior mean values for the fit.
        priorsigma: array_like, optional, default = None
            For constrained fits. Prior error bars for the fit.
        detailedInfo : bool, optional. default = False
            If True, return also the covariance matrix of the fit parameters and logGBF

        Returns
        -------
        params:
            The final fit parameters.
        params_err:
            The error of this fit parameters.
        chidof:
            chi^2/dof.
        stats:
            If detailedInfo, return dictionary with logGBF, BAIC, pcov.
        """

        if ('curve_fit' in algorithms) and self._errorAlg=='hessian':
            self._errorAlg='propagation'
            logger.warn('hessian strategy not yet compatible with curve_fit; switching strategy to propagation.')

        # Prior values are a good point for starting the fit, if you don't have another guess.
        if (priorval is not None) and (start_params is None):
            start_params = priorval

        # For a fit with only one parameter, we also accept a scalar. Check if this is the case.
        if start_params is not None:
            self._saved_params = envector(start_params)

        # Initialize prior values.
        checkPrior(priorval,priorsigma)
        if priorsigma is not None:
            self._priorsigma = np.array(priorsigma)
            self._priorval = np.array(priorval)

        # See if a standard function evaluation works
        try:
            _test = self.wrap_func(self._xdata, self._saved_params)
        except Exception as e:
            logger.TBError('Fit function must have signature func(xdata,params,args). Got exception:',e)
        checkEqualLengths(_test, self._xdata)

        resultSummary  = parallel_function_eval( self._tryAlgorithm, algorithms, nproc=self._nproc )
        all_params     = [row[0] for row in resultSummary]
        all_fit_errors = [row[1] for row in resultSummary]
        all_chi2       = [row[2] for row in resultSummary]
        all_pcov       = [row[3] for row in resultSummary]
        all_except     = [row[4] for row in resultSummary]

        # Check to see if all the fits failed
        if np.all(np.array(all_chi2) == np.inf): 
            for i, algorithm in enumerate(algorithms):
                logger.TBFail(algorithm+": ",all_except[i])
            logger.TBError("No algorithm converged. See above list of exceptions.")

        # Find the smallest chi^2
        min_ind = np.argmin(all_chi2)
        if len(algorithms) > 1:
            logger.details("Choosing", algorithms[min_ind], "with chi^2 =", all_chi2[min_ind])

        # Store results internally
        self._saved_params = all_params[min_ind]
        self._saved_pcov = all_pcov[min_ind]

        for i in range(len(algorithms)):
            logger.debug(algorithms[i]+':',all_params[i],all_fit_errors[i],all_chi2[i],all_pcov[i],all_except[i])
            
        dof = DOF(len(self._ydata),len(self._saved_params),self._priorsigma)    

        if dof <= 0:
            chidof = np.inf
        else:
            chidof = all_chi2[min_ind]/dof

        if detailedInfo:
            return ( np.copy(self._saved_params),
                     all_fit_errors[min_ind],
                     chidof,
                     {  'logGBF'   : logGBF(self._xdata, self._ydata, self._cov, self._func, self._args, self._saved_params,
                                           prior=self._priorval,priorsigma=self._priorsigma),
                        'pcov'     : np.copy(self._saved_pcov),
                        'chi2'     : all_chi2[min_ind],
                        'chi2data' : chisquare(self._xdata, self._ydata, self._cov, self._func, self._args, self._saved_params,
                                               None, None),
                        'BAIC'     : BAIC(self._xdata, self._ydata, self._cov, self._func, self._args, self._saved_params),
                        'AIC'      : AIC(self._xdata, self._ydata, self._cov, self._func, self._args, self._saved_params,
                                         self._priorval, self._priorsigma),
                        'AICc'     : AICc(self._xdata, self._ydata, self._cov, self._func, self._args, self._saved_params,
                                          self._priorval, self._priorsigma),
                        'DOF'      : dof,
                        'Q'        : goodnessOfFit(dof,all_chi2[min_ind])
                     }
                   )
        else:
            return ( np.copy(self._saved_params),
                     all_fit_errors[min_ind],
                     chidof
                   )


    def do_fit(self, algorithm="curve_fit", **kwargs):
        """ Same as try_fit but with only one algorithm. """
        return self.try_fit([algorithm], **kwargs)


    def save_func(self, filename, domain = None, no_error=False, header=None, npoints=1000, **kwargs):
        """ Save fit data to table. """
        domain = self._autoDomain(domain)
        params_err = self._saved_pcov
        params = self._saved_params
        def func(x, params):
            return self._func(x, params, *self._args)
        if no_error:
            save_func(func, filename, domain=domain, args=params, header=header, npoints=npoints,**kwargs)
        else:
            def grad(x, params):
                return np.asarray(self.grad(x, params))
            save_func(func, filename, domain=domain, args=params, args_err=params_err, grad=grad,
                      header=header, npoints=npoints, **kwargs)


    def plot_fit(self, domain = None, no_error = False, **kwargs):
        """ Plot the fit function. """
        logger.debug('Plotting fit.')
        domain = self._autoDomain(domain)
        if not no_error:
            plot_func(self._func, domain=domain, params=self._saved_params, params_err=self._saved_pcov, 
                      grad=self._grad, args=self._args, **kwargs)
        else:
            plot_func(self._func, domain=domain, params=self._saved_params, args=self._args, **kwargs)


    def plot_data(self, **kwargs):
        """ Plot the fit data. """
        logger.debug('Plotting fit data.')
        if self._cov is not None:
            sigma = np.sqrt( np.diag(self._cov) )
            plot_dots(self._xdata, self._ydata, sigma, **kwargs)
        else:
            plot_dots(self._xdata, self._ydata, **kwargs)


    def plot_cor(self, title = 'Data correlation matrix', xlabel='$x_j$', ylabel='$x_i$'):
        """ Plot the correlation matrix of the input data. """

        ncov = self._fit_cor

        if title is not None:
            plt.title(title)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        xrange = np.arange(len(ncov) + 1)
        yrange = np.arange(len(ncov) + 1)

        plt.pcolormesh(xrange, yrange, ncov, cmap="Blues")
        plt.gca().invert_yaxis()
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=12)


    def plot_eig(self, xlabel="$i$", ylabel="$\\lambda_i$", title = 'Eigenvalues of data correlation matrix'):

        vals, _ = np.linalg.eig(self._fit_cor)

        eig_real = np.real(vals)
        eig_imag = np.imag(vals)

        plot_bar(range(len(eig_real)), eig_real, color='#d32d11', label="real", alpha=0.7, title=title, xlabel=xlabel, ylabel=ylabel)
        if np.min(eig_imag) != 0:
            plot_bar(range(len(eig_imag)), eig_imag, color='#0081bf', label="imag", alpha=0.7, title=title, xlabel=xlabel, ylabel=ylabel)


def save_func(func, filename, domain, args=(), func_err=None, args_err=(), grad = None, header=None, 
              npoints=1000, **kwargs):

    checkDomain(domain)
    xmin = domain[0] 
    xmax = domain[1] 

    if xmin is None:
        for line in plt.gca().lines:
            xmin_new = np.min(line.get_xdata())
            if xmin is None:
                xmin = xmin_new
            if xmin_new < xmin:
                xmin = xmin_new
    if xmax is None:
        for line in plt.gca().lines:
            xmax_new = np.max(line.get_xdata())
            if xmax is None:
                xmax = xmax_new
            if xmax_new > xmax:
                xmax = xmax_new

    xdata = np.arange(xmin, xmax, (xmax - xmin) / npoints)
    ydata = func(xdata, args)

    createFilePath(filename)

    if func_err is not None:
        ydata_err = func_err(xdata, args_err)
        writeTable(filename,xdata,ydata,ydata_err,header=header,**kwargs)

    elif len(args_err) > 0:
        if grad is None:
            logger.warn("Used numerical derivative!")

        # Arguments that are part of the error propagation
        tmp_args = tuple(args)[0:len(args_err)]

        # Optional arguments that are constant and, therefore, not part of the error propagation
        tmp_opt = tuple(args)[len(args_err):]

        ydata_err = error_prop_func(xdata, func, tmp_args, args_err, grad=grad, args=tmp_opt)

        writeTable(filename,xdata,ydata,ydata_err,header=header,**kwargs)

    else:
        writeTable(filename,xdata,ydata,header=header,**kwargs)


def do_fit(func, xdata, ydata, edata=None, start_params=None, priorval=None, priorsigma=None,
           algorithm="curve_fit", detailedInfo=False, **kwargs):
    """ Wrapper to fitter initialization and the fit in one step. See above for arguments. """
    fit = Fitter(func, xdata, ydata, edata, **kwargs)
    return fit.do_fit(start_params=start_params, priorval=priorval, priorsigma=priorsigma, algorithm=algorithm,
                      detailedInfo=detailedInfo)


def try_fit(func, xdata, ydata, edata=None, start_params=None, priorval=None, priorsigma=None,
            algorithms=std_algs, detailedInfo=False, **kwargs):
    """ Wrapper to fitter initialization and the fit in one step. See above for arguments. For historical reasons
    algorithms has no default values here. """
    fit = Fitter(func, xdata, ydata, edata, **kwargs)
    return fit.try_fit(algorithms, start_params, priorval, priorsigma, detailedInfo=detailedInfo)


