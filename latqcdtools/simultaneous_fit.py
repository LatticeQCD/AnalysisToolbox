from latqcdtools.fitting import Fitter
import numpy as np


"""
Fit simultaneously different sets of data. For example in a continuum extrapolation.
Each set of data should be identified by an additional parameter (add_param).
A typical example would be Nt or the lattice spacing.
This parameters should be passed as the second argument to the fit function.

Parameters
----------
func:
    The function to be fitted
add_params:
    List of additional parameters that are needed to distinguish the data sets.
xdata:
    List of sets of xdata
ydata:
    List of sets of ydata
edata:
    List of sets of error data

Rest of parameters is analog to the standard fitting parameters in fitting.py
"""


def simultaneous_fit(func, add_params, xdata, ydata, edata = None, 
        func_sup_numpy = False, start_params = None, 
        args = (), expand = True, algorithms = ['levenberg'], ret_pcov = False, **kwargs):



    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)


    lengths = np.array([len(x) for x in xdata])
    xdata = np.concatenate(xdata)
    ydata = np.concatenate(ydata)
    if edata is None:
        edata = np.ones_like(ydata)
    else:
        edata = np.concatenate(edata)




    def fit_func(x, params):
        ret = np.zeros_like(ydata)
        for i, xset in enumerate(add_params):

            ind_dn = np.sum(lengths[:i])
            ind_up = ind_dn + lengths[i]

            if func_sup_numpy:
                if expand:
                    ret[ind_dn:ind_up] = func(x[ind_dn:ind_up], add_params[i],
                            *(tuple(params) + tuple(args)))
                else:
                    ret[ind_dn:ind_up] = func(x[ind_dn:ind_up], add_params[i], params, *args)
            else:
                for j in range(ind_dn, ind_up):
                    if expand:
                        ret[j] = func(x[j], add_params[i], *(tuple(params) + tuple(args)))
                    else:
                        ret[j] = func(x[j], add_params[i], params, *args)
        return ret


    fitter = Fitter(fit_func, xdata, ydata, edata, expand = False,
            func_sup_numpy = True, **kwargs)
    return fitter.try_fit(algorithms = algorithms, start_params = start_params,
            ret_pcov = ret_pcov)




