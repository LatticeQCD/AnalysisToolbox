# 
# optimize.py
# 
# H. Sandmeyer, D. Clarke
# 
# Some methods related to minimization and Levenberg fits.
#
import warnings, multiprocessing
import scipy.optimize as opt
from scipy.optimize import newton_krylov, fsolve, root
import numpy as np
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("error", category=RuntimeWarning)
from scipy.optimize.nonlin import NoConvergence
import latqcdtools.base.logger as logger
from latqcdtools.base.check import DivideByZeroError, InvalidValueError
from latqcdtools.base.utilities import unvector


# This is the base list of exceptions. If encountered, we treat the solve as unreliable.
opt_exceptions = (NoConvergence, FloatingPointError, ValueError, RuntimeWarning, DivideByZeroError, InvalidValueError)


def timeout(func, args=(), kwargs={}, timeout_duration=300):
    """ Exit if the fit is taking too long ."""
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    def wrap_func(*wrap_args):
        ret = func(*wrap_args[0], **wrap_args[1])
        wrap_args[2]["ret"] = ret

    p = multiprocessing.Process(target = wrap_func, args = (args, kwargs, return_dict))
    p.start()

    p.join(timeout_duration)

    # If thread is still active
    if p.is_alive():
        p.terminate()
        p.join()
        raise TimeoutError("Time out for " + str(func))
    return return_dict["ret"]


def persistentSolve(LHS, guess, tol=1e-8, maxiter=200):
    """ Attempt to solve LHS==0 using, in this order, SciPy's newton_krylov, fsolve, and root. """
    try:
        logger.debug("Trying newton_krylov.")
        solution = unvector( newton_krylov(LHS, guess, f_tol=tol, inner_maxiter=maxiter) )
    except opt_exceptions:
        try:
            logger.debug("Trying fsolve.")
            solution = unvector( fsolve(LHS, guess, xtol=tol, maxfev=maxiter) )
        except opt_exceptions:
            try:
                logger.debug("Trying root.")
                solution = unvector( root(LHS, guess, tol=tol) )
                if type(solution).__name__ == "OptimizeResult":
                    raise RuntimeWarning('RuntimeWarning: Bad progress on root solve.')
            except Exception as e:
                raise e
    return unvector(solution)


def minimize(func, jack=None, hess=None, start_params=None, tol=1e-12, maxiter=10000, algorithm=None):

    args = (func, start_params)

    kwargs = {'method': algorithm,
                 'tol': tol}

    if algorithm == "BFGS":
        kwargs['jac'] = jack
        kwargs['options'] = {'gtol': tol, 'maxiter': maxiter}
        logger.details('I am using a gradient.')

    elif algorithm in ["TNC","SLSQP","L-BFGS-B","CG"]:
        kwargs['jac'] = jack
        kwargs['options'] = {'maxiter': maxiter}
        logger.details('I am using a gradient.')

    elif algorithm in ["COBYLA","Nelder-Mead"]:
        kwargs['options'] = {'maxiter': maxiter}

    elif algorithm == "Powell":
        kwargs['options'] = {'xtol': tol, 'ftol': tol, 'maxfev': maxiter}

    else:
        kwargs['jac'] = jack
        kwargs['hess'] = hess
        kwargs['options'] = {'maxiter': maxiter}
        logger.details('I am using a gradient and hessian.')

    # At least COBYLA sometimes gets stuck in an endless loop.
    res = timeout(opt.minimize, args=args, kwargs=kwargs, timeout_duration=100)

    params = res.x
    nfev = res.nfev
    if not res.success:
        logger.details(algorithm, res.message)
        raise ValueError(algorithm + ": Minimization did not converge!")

    try:
        params[0]
    except:
        if isinstance(params, (np.ndarray, np.generic)):
            # Somehow numpy 0D arrays have to be converted to scalar explicitly.
            params = [np.asscalar(params)]
        else:
            params = [params]

    return params, nfev
