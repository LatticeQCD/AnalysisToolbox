

import numpy as np



def diff_deriv(x, func, args = (), expand = False, eps = None):
    if eps is None:
        eps = best_eps(x)
    up = x + eps
    down = x - eps
    return (func(up, *args) - func(down, *args)) / (2*eps)





# Hessian using difference quotient
def diff_hess(params, func, args = (), eps = None, expand = False):
    if expand:
        wrap_func = lambda params: func(*(tuple(params) + tuple(args)))
        return diff_hess(params, wrap_func, eps = eps, expand = False)

    # This has to be a list, as we might put in arrays, if params is higher dimensional
    ret = [ [0.0 for i in range(len(params))] for j in range(len(params))]
    # convert to float
    up = np.array(params, dtype = float)
    down = np.array(params, dtype = float)
    upup = np.array(params, dtype = float)
    updown = np.array(params, dtype = float)
    downdown = np.array(params, dtype = float)
    downup = np.array(params, dtype = float)

    for i in range(len(params)):
        if eps is None:
            epsi = best_eps(params[i])
        else:
            epsi = eps
        up[i] += epsi
        upup[i] = up[i]
        updown[i] = up[i]
        down[i] -= epsi
        downdown[i] = down[i]
        downup[i] = down[i]
        ret[i][i]=(func(up,*args)+func(down,*args)-2*func(params,*args))/(4*(epsi/2.0)**2)

        for j in range(i):
            if eps is None:
                epsj = best_eps(params[j])
            else:
                epsj=eps

            upup[j]+=epsj
            updown[j]-=epsj
            downdown[j]-=epsj
            downup[j]+=epsj

            ret[i][j]=(func(upup,*args) + func(downdown,*args)
                    -func(updown,*args) - func(downup,*args)) / (4*epsi*epsj)
            ret[j][i]=ret[i][j]

            upup[j]=params[j]
            updown[j]=params[j]
            downdown[j]=params[j]
            downup[j]=params[j]

        up[i] = params[i]
        upup[i] = up[i]
        updown[i] = up[i]
        down[i] = params[i]
        downdown[i] = down[i]
        downup[i] = down[i]


    return np.array(ret)





def diff_fit_grad(x, params, func, args = (), eps = None, expand = False):
# For fitting or plotting we expect the first argument of func to be x instead of params.
# Therefore we have to change the order using this wrapper
    if expand:
        f = lambda p: func(x, *(tuple(p) + tuple(args)))
    else:
        f = lambda p: func(x, p, *args)
    return diff_grad(params, f, eps = eps, expand = False)




def diff_fit_hess(x, params, func, args = (), eps = None, expand = False):
# For fitting or plotting we expect the first argument of func to be x instead of params.
# Therefore we have to change the order using this wrapper
    if expand:
        f = lambda p: func(x, *(tuple(p) + tuple(args)))
    else:
        f = lambda p: func(x, p, *args)
    return diff_hess(params, f, eps = eps, expand = False)





# Implementation for usage with algopy which can handle numpy operations but is very slow and buggy




# alg_hess using algopy
def alg_hess(x, func, args = ()):
    if not algopy_avail:
        raise ImportError("Algopy not installed. Please install algopy or use diff_hess")
    return alg_jac(x, alg_grad, args = (func, args))


def alg_fit_grad(x, params, func, args = (), expand = False):
# For fitting or plotting we expect the first argument of func to be x instead of params.
# Therefore we have to change the order using this wrapper
    if expand:
        f = lambda p: func(x, *(tuple(p) + tuple(args)))
    else:
        f = lambda p: func(x, p, *args)
    return alg_grad(params, f)

def alg_fit_hess(x, params, func, args = (), expand = False):
# For fitting or plotting we expect the first argument of func to be x instead of params.
# Therefore we have to change the order using this wrapper
    if expand:
        f = lambda p: func(x, *(tuple(p) + tuple(args)))
    else:
        f = lambda p: func(x, p, *args)
    return alg_hess(params, f)


