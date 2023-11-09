# 
# optimize.py
# 
# H. Sandmeyer, D. Clarke
# 
# Some methods related to minimization and Levenberg fits.
#

import scipy.optimize as opt
from scipy.optimize import newton_krylov, fsolve, root
from scipy.optimize.nonlin import NoConvergence
import latqcdtools.base.logger as logger
from latqcdtools.base.check import DivideByZeroError, InvalidValueError, checkType
from latqcdtools.base.utilities import envector 


# This is the base list of exceptions. If encountered, we treat the solve as unreliable.
opt_exceptions = (NoConvergence, FloatingPointError, ValueError, RuntimeWarning, DivideByZeroError, InvalidValueError)

# This order is not necessarily optimal.
root_methods = ['hybr','lm','broyden1','broyden2','anderson','diagbroyden','krylov']


def persistentSolve(LHS, guess, tol=1e-8, maxiter=300):
    """ Attempt to solve LHS==0 using, in this order, SciPy's newton_krylov, fsolve, and root. """
    checkType(tol,float)
    checkType(maxiter,int)
    try:
        logger.debug("Trying newton_krylov.")
        return newton_krylov(LHS, guess, f_tol=tol, inner_maxiter=maxiter)
    except opt_exceptions:
        try:
            logger.debug("Trying fsolve.")
            return fsolve(LHS, guess, xtol=tol, maxfev=maxiter)
        except opt_exceptions:
            try:
                logger.debug("Trying root.")
                for method in root_methods:
                    try:  
                        return root(LHS, guess, tol=tol, method=method).x
                    except:
                        continue
            except:
                raise NoConvergence


def minimize(func, jack=None, hess=None, start_params=None, tol=1e-12, maxiter=10000, algorithm=None):

    kwargs = {'method': algorithm, 'tol': tol}

    logger.details('Trying',algorithm,'with maxiter=',maxiter)

    if algorithm in ["TNC","SLSQP","L-BFGS-B","CG"]:
        kwargs['jac'] = jack
        kwargs['options'] = {'maxiter': maxiter}

    elif algorithm in ["COBYLA","Nelder-Mead"]:
        kwargs['options'] = {'maxiter': maxiter}

    elif algorithm == "Powell":
        kwargs['options'] = {'xtol': tol, 'ftol': tol, 'maxfev': maxiter}

    else:
        kwargs['jac'] = jack
        kwargs['hess'] = hess
        kwargs['options'] = {'maxiter': maxiter}

    try:
        res = opt.minimize(func, start_params, **kwargs)
    except Exception as e:
        raise e

    logger.debug(res)

    params = res.x
    if not res.success:
        logger.details(algorithm, res.message)
        raise ValueError('Minimizer failed.')

    params = envector(params)

    return params
