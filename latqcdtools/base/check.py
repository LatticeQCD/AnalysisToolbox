#
# check.py
#
# H. Sandmeyer, D. Clarke
#
# A collection of methods for comparing results, formatting test output, and handling numpy exceptions.
#

import math
import numpy as np
import latqcdtools.base.logger as logger


class DivideByZeroError(Exception): pass
class UnderflowError(Exception): pass
class InvalidValueError(Exception): pass


def err_handler(err, flag):
    if flag == 1:
        raise DivideByZeroError(err)
    elif flag == 2:
        raise OverflowError
    elif flag == 4:
        raise UnderflowError(err)
    elif flag == 8:
        raise InvalidValueError(err)


np.seterrcall(err_handler)
np.seterr(all='call')


def rel_check(a, b, prec = 1e-6, abs_prec = 1e-14):
    """ Check whether two values are equal. Use especially when comparing to 0. """
    return math.isclose(a, b, rel_tol = prec, abs_tol = abs_prec)


def rel_checkArrayScalar(arr, scal, prec = 1e-6, abs_prec = 1e-14):
    """ Return a boolean array that checks element-wise whether a numpy array arr is equal to a scalar scal. """
    comparisonArray = scal * np.ones(arr.shape)
    return np.isclose(arr, comparisonArray, prec, abs_prec)


def print_results(res, res_true, res_err = None, res_err_true = None, text = "", prec = 1e-10):
    """ Compares element-by-element the results of res with res_true. (Does the same with res_err and res_err_true,
    if you like.) Carries out with precision prec. """
    test = True

    try:
        res[0]
    except (IndexError, TypeError):
        res = [res]

    try:
        res_true[0]
    except (IndexError, TypeError):
        res_true = [res_true]

    if res_err is not None:
        try:
            res_err[0]
        except (IndexError, TypeError):
            res_err = [res_err]

    if res_err_true is not None:
        try:
            res_err_true[0]
        except (IndexError, TypeError):
            res_err_true = [res_err_true]

    for i in range(len(res)):
        if not rel_check(res[i], res_true[i], prec):
            test = False
            print("res[" + str(i) + "] = " + str(res[i])
                  + " != res_true[" + str(i) + "] = " + str(res_true[i]))
        if res_err is not None and res_err_true is not None:
            if not rel_check(res_err[i], res_err_true[i], prec):
                test = False
                print("res_err[" + str(i) + "] = " + str(res_err[i])
                      + " != res_err_true[" + str(i) + "] = " + str(res_err_true[i]))

    if test:
        logger.TBPass(text)
    else:
        logger.TBFail(text)