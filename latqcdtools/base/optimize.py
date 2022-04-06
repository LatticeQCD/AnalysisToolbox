# 
# optimize.py
# 
# H. Sandmeyer
# 
# Some methods related to minimization and Levenberg fits.
#
import multiprocessing
import scipy.optimize as opt
import numpy as np
from numpy.linalg import inv
from latqcdtools.base.num_deriv import diff_grad, diff_hess
import latqcdtools.base.logger as logger


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
        # Terminate
        p.terminate()
        p.join()
        raise TimeoutError("Time out for " + str(func))
    return return_dict["ret"]


def add_diag(mat, value):
    mat += value * np.diag(np.diag(mat))
    return mat


# -------------------------------------------------------------------------------------------------------- LEVENBERG FIT


def levenberg_step(func, grad_func, hess_func, param, lamb, args=()):
    if grad_func is None:
        grad = diff_grad(param, func, args)
    else:
        grad = grad_func(param, *args)

    if hess_func is None:
        hess = diff_hess(param, func, args)
    else:
        hess = hess_func(param, *args)
    tmp = add_diag(hess, lamb)
    param_new = param - inv(tmp).dot(grad)
    old_func_val = func(param, *args)
    diff = func(param_new, *args) - old_func_val
    if diff < 0:
        lamb /= 10
        return param_new, lamb, abs(diff) / old_func_val
    else:
        lamb *= 10
        return param, lamb, abs(diff) / old_func_val


def levenberg(param, func, grad_func=None, hess_func=None, eps=1e-12, max_itt=10000, print_output=False, args=()):
    diff = 1
    lamb = 1e-6
    i = 0
    while diff > eps:
        param, lamb, diff = levenberg_step(func, grad_func, hess_func, param, lamb, args = args)
        if np.isnan(diff):
            raise ValueError("Levenberg: Invalid value during function evaluation")
        if print_output:
            print("Iteration:", i, "Current parameter:", param, "Lambda:", lamb, "Relative difference", diff)
        i += 1
        if i > max_itt:
            raise ValueError("Levenberg: Number of iterations exceeded")
    return param, i


# --------------------------------------------------------------------------------------------------------- MINIMIZATION


def minimize(func, jack=None, hess=None, start_params=None, tol=1e-12, maxiter=10000, algorithm=None):
    if algorithm == "levenberg":
        args = (start_params, func, jack, hess)
        kwargs = {'eps': tol, 'max_itt': maxiter}
        params, nfev = levenberg(*args, **kwargs)

    else:
        args = (func, start_params)

        if algorithm == "BFGS":
            kwargs = {'method': algorithm,
                      'jac': jack,
                      'tol': tol,
                      'options': {'gtol': tol, 'maxiter': maxiter}}

        elif algorithm == "TNC":
            kwargs = {'method': algorithm,
                      'jac': jack,
                      'tol': tol,
                      'options': {'maxiter': maxiter}}

        elif algorithm == "COBYLA":
            kwargs = {'method': algorithm,
                      'tol': tol,
                      'options': {'maxiter': maxiter}}

        elif algorithm == "SLSQP":
            kwargs = {'method': algorithm,
                      'jac': jack,
                      'tol': tol,
                      'options': {'maxiter': maxiter}}

        elif algorithm == "L-BFGS-B":
            kwargs = {'method': algorithm,
                      'jac': jack,
                      'tol': tol,
                      'options': {'maxiter': maxiter}}

        elif algorithm == "Powell":
            kwargs = {'method': algorithm,
                      'tol': tol,
                      'options': {'xtol': tol, 'ftol': tol,
                                  'maxfev': maxiter}}

        elif algorithm == "Nelder-Mead":
            kwargs = {'method': algorithm,
                      'tol': tol,
                      'options': {'maxiter': maxiter}}

        else:
            kwargs = {'method': algorithm,
                      'jac': jack,
                      'hess': hess,
                      'tol': tol,
                      'options': {'maxiter': maxiter}}

        # At least COBYLA sometimes gets stuck in an endless loop.
        res = timeout(opt.minimize, args=args, kwargs=kwargs, timeout_duration=100)

    if algorithm != "levenberg":
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
