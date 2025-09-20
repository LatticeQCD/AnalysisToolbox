# 
# continuumExtrap.py                                                               
# 
# D. Clarke
# 
# A simple module for performing extrapolations to fit functions.
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.plotting import plt
from latqcdtools.statistics.fitting import Fitter
from latqcdtools.base.speedify import DEFAULTTHREADS
from latqcdtools.base.check import checkType


def _powerSeries(x,coeffs):
    """ 
    The default fit form for a continuum extrapolation is a power series in a^2.

    Args:
        x (array-like): a^2 or 1/Nt^2 
        coeffs (array-like): power series coefficients 

    Returns:
        array-like: power series in x 
    """
    result = 0.
    for i in range(len(coeffs)):
        result += coeffs[i]*x**i
    return result


class Extrapolator(Fitter):


    def __init__(self, x, obs, obs_err, order=1, xtype="a", error_strat='propagation', ansatz=None, nproc=DEFAULTTHREADS,
                 tol=1e-12, max_fev=None):
        """ 
        A framework for doing continuum limit extrapolations.

        Args:
            x (array-like): a data or Nt data. (It will get squared in here.) 
            obs (array-like)
            obs_err (array-like)
            ansatz (func, optional): continuum-limit fit ansatz. Power series in a^2 by default.
            order (int, optional): order of the power series. Defaults to 1.
            xtype (str, optional): choose to fit a data or Nt data. Defaults to "a".
            error_strat (str, optional): calculate errors using error propagation or augmented chi^2. Defaults to 'propagation'.
            nproc (int, optional): number of processors for fitting. Defaults to DEFAULTTHREADS.
            tol (float, optional): tolerance for the minimization. Defaults to 1e-12
            max_fev (int, optional): maximum number of iterations. Defaults to 10000
        """
        checkType('int',order=order)
        checkType(np.ndarray,x=x)
        checkType(np.ndarray,obs=obs)
        checkType(np.ndarray,obs_err=obs_err)

        self._order              = order
        self._triedExtrapolation = False
        self._ansatz             = ansatz

        if xtype == "a":
            x = x**2
        elif xtype == "Nt":
            x = 1/x**2
        else:
            logger.TBRaise('Unknown xtype',xtype)
        if order<1:
            logger.TBRaise('Please input order > 1.')

        if ansatz is None:
            ansatz = _powerSeries
        else:
            if self._order != 1:
                logger.warn('Not using a power series ansatz, but still using custom order.')
            logger.debug('received ansatz',ansatz)

        Fitter.__init__(self, ansatz, x, obs, obs_err, norm_err_chi2=False, error_strat=error_strat, nproc=nproc,
                        tol=tol,max_fev=max_fev)


    def __repr__(self) -> str:
        return "Extrapolator"


    def extrapolate(self,start_coeffs=None,prior=None,priorsigma=None,detailedInfo=False,show_results=False):
        """ 
        Carry out the extrapolation.

        Args:
            start_coeffs (array-like, optional): your guess for starting parameters. Defaults to None.
            prior (array-like, optional): Bayesian priors. Defaults to None.
            priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
            detailedInfo (bool, optional): Do you want information like AIC? Defaults to False.
            show_results (bool, optional): Print fit results to screen? Defaults to False.
            
        Returns:
            (array-like, array-like, float, float): fit result, its error, chi^2/d.o.f., and logGBF if relevant.
        """
        if start_coeffs is None:
            if self._ansatz is not None:
                logger.TBRaise('You need to provide start_coeffs if you use an ansatz.')
            coeffs = np.ones(self._order+1) 
        else:
            checkType(np.ndarray,start_coeffs=start_coeffs)
            coeffs=start_coeffs
        self._triedExtrapolation = True
        if prior is None:
            self._result, self._result_err, self._chidof, self._stats = self.try_fit(start_params=coeffs, algorithms=["curve_fit"], 
                                                                                     detailedInfo=True, show_results=show_results)
        else:
            self._result, self._result_err, self._chidof, self._stats = self.try_fit(start_params=coeffs, priorval=prior, 
                                                                                     priorsigma=priorsigma, show_results=show_results, 
                                                                                     algorithms=["Nelder-Mead"], detailedInfo=True)
        if detailedInfo:
            return self._result, self._result_err, self._chidof, self._stats
        else:
            return self._result, self._result_err, self._chidof


    def plot(self,**kwargs):
        """ 
        Add extrapolation to plot. Accepts the same kwargs as Fitter.plot_fit. 
        """
        if not self._triedExtrapolation:
            logger.TBRaise("Can't plot an extrapolation without having extrapolated first...")
        domain=(1e-8,np.max(self._xdata))
        self.plot_data(**kwargs)
        kwargs['label']=None
        self.plot_fit(domain,**kwargs)


    def save(self,filename,header):
        if not self._triedExtrapolation:
            logger.TBRaise("Can't save an extrapolation without having extrapolated first...")
        domain=(1e-8,np.max(self._xdata))
        self.save_func(filename,domain,header=header)


def continuumExtrapolate(x,obs,obs_err,order=1,show_results=False,plot_results=False,prior=None, start_coeffs=None,priorsigma=None,
                         error_strat='propagation',xtype="a",nproc=DEFAULTTHREADS,detailedInfo=False):
    """ 
    A convenience wrapper for the Extrapolator. 
    """
    ext = Extrapolator(x, obs, obs_err, xtype=xtype, order=order, error_strat=error_strat, nproc=nproc)
    result = ext.extrapolate(start_coeffs=start_coeffs, prior=prior, priorsigma=priorsigma, detailedInfo=detailedInfo, show_results=show_results)
    if plot_results:
        if xtype=="a":
            xlabel = "$a^2$"
        else:
            xlabel = "$1/N_\\tau^2$"
        ext.plot(xlabel=xlabel)
        plt.show()
    return result
