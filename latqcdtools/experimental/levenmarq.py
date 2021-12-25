import numpy as np
from numpy.linalg import inv
from latqcdtools.num_deriv import alg_grad, diff_grad, alg_hess, diff_hess


def func(a):
    return a[0]**2 + a[1]**2


def grad_func(a):
    return [2 * a[0], 2 * a[1]]


def hess_func(a):
    return np.array([[2., 0.], [0., 2.]])


def add_diag(mat, value):
    mat += value * np.diag(np.diag(mat))
    return mat


def levenberg_step(func, grad_func, hess_func, param, lamb, use_alg, args=()):
    if grad_func is None:
        if use_alg:
            grad = alg_grad(param, func, args)
        else:
            grad = diff_grad(param, func, args)
    else:
        grad = grad_func(param, *args)

    if hess_func is None:
        if use_alg:
            hess = alg_hess(param, func, args)
        else:
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


def levenberg(param, func, grad_func=None, hess_func=None, eps=1e-12, max_itt=10000,
        print_output=False, use_alg = False, args=()):
    diff = 1
    lamb = 1e-6
    i = 0

    while diff > eps:
        param, lamb, diff = levenberg_step(func, grad_func, hess_func,
                param, lamb, use_alg, args = args)
        if np.isnan(diff):
            raise ValueError("Levenberg: Invalid value during function evaluation")

        if print_output:
            print("Iteration:", i, "Current parameter:", param, "Lambda:", lamb,
                    "Relative difference", diff)
        i += 1
        if i > max_itt:
            raise ValueError("Levenberg: Number of iterations exceeded")
    return param, i


def test_levenberg():
    print(levenberg(func, grad_func, hess_func, [10, 10]))

