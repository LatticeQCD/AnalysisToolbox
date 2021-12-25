import numpy as np
from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline
import latqcdtools.logger as logger
import ctypes, os
from scipy import interpolate
from latqcdtools.simultaneous_fit import simultaneous_fit
from latqcdtools.spline import constraint_spline, get_num_spline_coeffs, spline_ind_sel


def spline(x, xdata, ydata, npoints = 1000):
    xmin=min(xdata)
    xmax=max(xdata)
    if x is None:
        x = np.linspace(xmin, xmax, npoints)
    tck = interpolate.splrep(xdata, ydata, s=0)
    ynew = interpolate.splev(x, tck, der=0)
    return ynew


def bspline_sci(xdata, ydata, ydata_err=None, xspline=None, npoints=1000, order = 3, s = None):
    spl = UnivariateSpline(xdata, ydata, ydata_err, k = order, s = s)
    xmin = min(xdata)
    xmax = max(xdata)
    if xspline is None:
        xspline=np.linspace(xmin, xmax, npoints)
    return xspline, spl(xspline)


def LSQbspline_sci(xdata, ydata, ydata_err = None, knots=None, ncoeffs = 10, order = 3,
        xspline = None, npoints = 1000):
    if knots is None:
        knots = np.linspace(np.max(xdata), np.min(xdata), ncoeffs)
    spl = LSQUnivariateSpline(xdata, ydata, knots, ydata_err, k = order)
    xmin = min(xdata)
    xmax = max(xdata)
    if xspline is None:
        xspline=np.linspace(xmin, xmax, npoints)
    return xspline, spl(xspline)


""" Bspline using gsl c library

INPUT:

    xdata
    ydata
    edata:   error of ydata
    xmin:    xmin of data used for spline calculation
    xmax:    xmax of data used for spline calculation
    sp_xmin: xmin for output data
    sp_xmax: xmax for output data
    xspline: instead of linspace between sp_xmin and sp_xmax use these x values for spline output
    ncoeffs: number of coeffecients. If None best choice is tuned automatically
    order:   order of the spline
    npoints: number of output points if not provided xspline
    s:       smoothing factor for automatic choice of ncoeffs.

OUTPUT:

    Coordinates (xout,yout,youterr) that make the graph of the spline.

"""
def bspline(xdata, ydata, edata = None, xmin = None, xmax = None, sp_xmin = None,
        xspline = None, sp_xmax = None, ncoeffs = None, order = 3,
        npoints = 1000, s = 1.1, silent=False, maxCoeffPerc=0.7):
    if xmin is None:
        xmin = np.min(xdata)
    if xmax is None:
        xmax = np.max(xdata)

    shared_lib_path = os.path.dirname(os.path.realpath(__file__)) + "/../bspline_build/libbspline_py.so"

    if edata is None:
        edata = np.ones_like(xdata)

    if xspline is None:
        xout = np.linspace(xmin, xmax, npoints)
    else:
        xout = xspline

    xout = np.sort(np.array(xout, dtype=float))
    yout = np.copy(xout)
    yout_err = np.copy(xout)

    xdata = np.array(xdata, dtype=float)
    ydata = np.array(ydata, dtype=float)
    edata = np.array(edata, dtype=float)

    ind = (xdata <= xmax) & (xdata >= xmin)

    xdata = xdata[ind]
    ydata = ydata[ind]
    edata = edata[ind]

    try:
        lib = ctypes.cdll.LoadLibrary(shared_lib_path)
    except OSError:
        print('You must compile bspline before usage. Go to ../bspline_build,'
                'type "cmake ../bspline; make"')
        quit(-1)
    bspline_interpolate = lib.bspline_interpolate

    if ncoeffs is None:
        diffs = []
        inds = []
        for nc in range(order + 1, int(len(xdata)*maxCoeffPerc)):
            res_x, res_y, res_e = bspline(xdata, ydata, edata,
                    ncoeffs = nc, order = order, xspline = xdata)
            chi_dof = np.sum((res_y - ydata)**2 / edata**2) / (len(xdata) - order)

            if chi_dof < s:
                ncoeffs = nc
                if not silent:
                    print("Bspline: Choose ncoeff = ", nc, "with chi/d.o.f =", chi_dof)
                break

        if ncoeffs is None:
            ncoeffs = int(len(xdata)*maxCoeffPerc)
            if not silent:
                print("Bspline: Could not fullfill chi/d.o.f <", s, ".")
                print("Choose ncoeff = length(xdata) = ", ncoeffs, "with chi/d.o.f =", chi_dof)

    bspline_interpolate(ctypes.c_void_p(xdata.ctypes.data), ctypes.c_void_p(ydata.ctypes.data),
                        ctypes.c_void_p(edata.ctypes.data), ctypes.c_int(len(xdata)), 
                        ctypes.c_void_p(xout.ctypes.data), ctypes.c_void_p(yout.ctypes.data), 
                        ctypes.c_void_p(yout_err.ctypes.data), ctypes.c_int(len(xout)), 
                        ctypes.c_int(order), ctypes.c_int(ncoeffs))
    return xout, yout, yout_err


"""
Calls a simultaneous_fit using constraint splines as fitting function.
"""
def multi_spline_fit(fit_func, nparams, add_params, xdata, ydata, edata = None, knots = None,
        nknots = 4, ncoeffs = None, order = 3, start_params = None, constraints = None,
        nsteady_deriv = None, randomization_factor = 0,
        fit_knots = False, base_point = 0, **kwargs):

    if constraints is None:
        constraints = []

    if nsteady_deriv is None:
        nsteady_deriv = order - 1

    if knots is not None:
        nknots = len(knots)

    if start_params is not None:
        # Flatten start_params if not already done
        try:
            start_params[0][0]
            start_params = np.concatenate(start_params)
        except (IndexError, TypeError):
            pass

    xmin = np.inf
    xmax = -np.inf

    for x in xdata:
        if np.min(x) < xmin:
            xmin = np.min(x)
        if np.max(x) > xmax:
            xmax = np.max(x)

    # Automatically distribute knots according to the density of points
    if knots is None:

        if randomization_factor > 0:
            knots = random_knots(xdata, nknots, constraints, randomization_factor)
        else:
            knots = auto_knots(xdata, nknots)

    knots = np.asarray(knots)

    # In many cases, constraint fits do not converge. Therefore we first perform
    # a non-constraint fit to estimate the start parameters.
    if start_params is None and (len(constraints) > 0):
        logger.progress("multi_spline_fit: No startparameters available,"
                "First unconstrained fit to estimate parameters")
        start_params = multi_spline_fit(fit_func, nparams, add_params, xdata, ydata,
                edata = edata, knots = knots,
                nknots = nknots, ncoeffs = ncoeffs, order = order,
                constraints = None,
                nsteady_deriv = nsteady_deriv, randomization_factor = 0,
                fit_knots = fit_knots, base_point = base_point, **kwargs)[1]

        remove_inds = spline_ind_sel(None, constraints, knots, order, nsteady_deriv)
        start_params = np.delete(start_params, remove_inds, axis = 0)
        start_params = np.concatenate(start_params)

    if ncoeffs is None:
        ncoeffs = get_num_spline_coeffs(order, nsteady_deriv, nknots, constraints)

    def func(x, add_param, params):
        coeffs = [fit_func(add_param,
            *params[nparams*j:nparams*j+nparams]) for j in range(ncoeffs)]

        # Also fit the knot positions
        if fit_knots:
            new_knots = [fit_func(add_param,
                *params[nparams*ncoeffs + nparams*j:nparams*ncoeffs + nparams*(j+1)])
                for j in range(nknots)]
            return constraint_spline(x, new_knots, coeffs, order, nsteady_deriv, nderiv = 0,
                    constraints = constraints, base_point = base_point)
        else:
            return constraint_spline(x, knots, coeffs, order, nsteady_deriv, nderiv = 0,
                    constraints = constraints, base_point = base_point)

    res, res_err, chi_dof = simultaneous_fit(func, add_params, xdata, ydata, edata,
            func_sup_numpy = True, expand = False, start_params = start_params, **kwargs)

    if fit_knots:
        # We can not return the knot positions if they are part of the fit parameters
        return (res.reshape(ncoeffs + fit_knots * nknots, nparams),
                res_err.reshape(ncoeffs + fit_knots * nknots, nparams), chi_dof)
    else:
        return (knots, res.reshape(ncoeffs + fit_knots * nknots, nparams),
                res_err.reshape(ncoeffs + fit_knots * nknots, nparams), chi_dof)


def auto_knots(xdata, nknots):
    try:
        flat_xdata = np.sort(np.concatenate(xdata))
    except ValueError:
        flat_xdata = np.sort(np.asarray(xdata))

    # to ensure no knot sits at a data position
    flat_xdata = np.unique(flat_xdata)

    jump_step = (len(flat_xdata) - 1) / (nknots + 1)

    knots = []
    for i in range(1, nknots + 1):
        x_lower = flat_xdata[int(np.floor(i * jump_step))]
        x_upper = flat_xdata[int(np.ceil(i * jump_step))]
        knots.append((x_lower + x_upper) / 2)

    return knots


def random_knots(xdata, nknots, constraints, randomization_factor = 1):
    try:
        flat_xdata = np.sort(np.concatenate(xdata))
    except ValueError:
        flat_xdata = np.asarray(xdata)

    
    sample_xdata = np.random.choice(flat_xdata,
            int(nknots + 1 + (1 - randomization_factor) * (len(flat_xdata) - nknots)),
            replace = False)

    # Retry if to many data points are removed by np.unique
    if len(np.unique(sample_xdata)) < nknots + 1:
        return random_knots(xdata, nknots, constraints, randomization_factor)


    ret = auto_knots(sample_xdata, nknots)

    # If per chance a knot is generated at a constraint position, we have to run again
    constraints = np.asarray(constraints)
    for constraint in constraints:
        if any(constraint[0] == ret):
            return random_knots(xdata, nknots, constraints, randomization_factor)
            
    return ret