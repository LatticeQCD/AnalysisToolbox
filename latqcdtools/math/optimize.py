# 
# optimize.py
# 
# H. Sandmeyer, D. Clarke
# 
# Some methods related to minimization and Levenberg fits.
#

import scipy.optimize as opt
from scipy.optimize import newton_krylov, fsolve, root
from scipy.optimize._nonlin import NoConvergence
import latqcdtools.base.logger as logger
from latqcdtools.base.check import DivideByZeroError, InvalidValueError, checkType
from latqcdtools.base.utilities import envector 


# This is the base list of exceptions. If encountered, we treat the solve as unreliable.
opt_exceptions = (NoConvergence, FloatingPointError, ValueError, RuntimeWarning, DivideByZeroError, InvalidValueError)


# These orders are not necessarily optimal.
_root_methods  = ['krylov','hybr','broyden1','broyden2','anderson','diagbroyden']
_solve_methods = ['newton_krylov','fsolve','root']


def solve(LHS,guess,tol=1e-8,maxiter=300,method='newton_krylov'):
    """ Wrapper for various methods to solve LHS==0. This is to simplify the interface.

    Args:
        LHS (func)
        guess: initial guess for solution 
        tol (real, optional): Solve tolerance. Defaults to 1e-8.
        maxiter (int, optional): Maximum iterations. Defaults to 300.
        method (str, optional): Defaults to 'newton_krylov'.

    """
    checkType(tol,"real")
    checkType(maxiter,int)
    if method=='newton_krylov':
        return newton_krylov(LHS, guess, f_tol=tol, maxiter=maxiter)
    elif method=='fsolve':
        return fsolve(LHS, guess, xtol=tol, maxfev=maxiter)
    elif method=='root':
        for root_method in _root_methods:
            try:  
                logger.debug("Trying root:",root_method)
                if root_method in ['hybr']:
                    res = root(LHS, guess, tol=tol, method=root_method, options={'maxfev':maxiter})
                else:
                    res = root(LHS, guess, tol=tol, method=root_method, options={'maxiter':maxiter})
                if res.success == True:
                    return res.x
                else:
                    raise NoConvergence(root_method+' did not converge') 
            except Exception as e:
                logger.debug("Hit exception:",str(e)+"; Terminating.")
                continue
        raise NoConvergence('No root method converged.')
    else:
        logger.TBError('Unrecognized method',method)


def persistentSolve(LHS, guess, tol=1e-8, maxiter=300):
    """ Attempt to solve LHS==0 using, in this order, SciPy's newton_krylov, fsolve, and root.

        Args:
        LHS (func)
        guess: initial guess for solution 
        tol (real, optional): Solve tolerance. Defaults to 1e-8.
        maxiter (int, optional): Maximum iterations. Defaults to 300.
    """
    checkType(tol,"real")
    checkType(maxiter,int)
    for method in _solve_methods:
        try:
            return solve(LHS,guess,tol=tol,maxiter=maxiter,method=method)
        except opt_exceptions as e:
            logger.debug("Hit exception:",e)
        except Exception as e:
            raise e
    raise NoConvergence('No method converged!') 


def minimize(func, jac=None, hess=None, start_params=None, tol=1e-12, maxiter=10000, algorithm=None):
    """ Wrapper for scipy.optimize.minimize. Helps that all algorithms have common syntax.

    Args:
        func: to-be-minimized function 
        jac (func, optional): Explicit Jacobian. Defaults to None.
        hess (func, optional): Explicit Hessian. Defaults to None.
        start_params (array-like, optional): Starting guess. Defaults to None.
        tol (float, optional): Solve tolerance. Defaults to 1e-12.
        maxiter (int, optional): Maximum number of solve iterations. Defaults to 10000.
        algorithm (str, optional): Solve algorithm. Defaults to Scipy default, usually BFGS.

    Raises:
        e: scipy.optimize exception 
        ValueError: minimizer fails to converge

    Returns:
        array-like: solution vector 
    """

    kwargs = {'method': algorithm, 'tol': tol}

    logger.details('Trying',algorithm,'with maxiter=',maxiter)

    if algorithm in ["SLSQP","L-BFGS-B","CG"]:
        kwargs['jac'] = jac
        kwargs['options'] = {'maxiter': maxiter}

    elif algorithm == "TNC":
        kwargs['jac'] = jac
        kwargs['options'] = {'maxfun': maxiter}

    elif algorithm in ["COBYLA","Nelder-Mead"]:
        kwargs['options'] = {'maxiter': maxiter}

    elif algorithm == "Powell":
        kwargs['options'] = {'xtol': tol, 'ftol': tol, 'maxfev': maxiter}

    else:
        kwargs['jac'] = jac
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

    return envector(params)
