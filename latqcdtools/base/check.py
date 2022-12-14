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
from latqcdtools.base.utilities import envector


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
    elif flag == 9:
        raise DivideByZeroError(err)
    else:
        logger.TBError('Encountered unknown exception',err,'with flag',flag)


np.seterrcall(err_handler)
np.seterr(all='call')


def checkType(instance, expectedType):
    """ Only useful if you expect one particular type. """
    if not isinstance(instance,expectedType):
        logger.TBError('Expected type',expectedType,'but received',type(instance),frame=3)


def rel_check(a, b, prec = 1e-6, abs_prec = 1e-14):
    """ Check whether two values are equal. Use especially when comparing to 0. """
    try:
        return math.isclose(a, b, rel_tol = prec, abs_tol = abs_prec)
    except TypeError:
        logger.TBError('Expected real numbers. Received a, b types = ',type(a),',',type(b))


def rel_checkArrayScalar(arr, scal, prec = 1e-6, abs_prec = 1e-14):
    """ Return a boolean array that checks element-wise whether a numpy array arr is equal to a scalar scal. """
    comparisonArray = scal * np.ones(arr.shape)
    return np.isclose(arr, comparisonArray, prec, abs_prec)


def print_results(res, res_true, res_err = None, res_err_true = None, text = "", prec = 1e-10, abs_prec = None):
    """ Compares element-by-element the results of res with res_true. (Does the same with res_err and res_err_true,
        if you like.) Carries out with precision prec. Use abs_prec for comparisons with zero. """
    test = True
    compareZero = False

    if abs_prec is not None:
        compareZero = True
    else:
        abs_prec = 1e-14

    res = envector(res)
    res_true = envector(res_true)
    if res_err is not None:
        res_err = envector(res_err)
    if res_err_true is not None:
        res_err_true = envector(res_err_true)

    for i in range(len(res)):
        if not rel_check(res[i], res_true[i], prec, abs_prec):
            test = False
            print("res[" + str(i) + "] = " + str(res[i]) + " != res_true[" + str(i) + "] = " + str(res_true[i]))
        if res_err is not None and res_err_true is not None:
            if not rel_check(res_err[i], res_err_true[i], prec, abs_prec):
                test = False
                print("res_err[" + str(i) + "] = " + str(res_err[i]) + " != res_err_true[" + str(i) + "] = " + str(res_err_true[i]))

    if test:
        if compareZero:
            logger.TBPass(text,'(abs_prec = %.2e)' % abs_prec)
        else:
            logger.TBPass(text,'(prec = %.2e)' % prec)
    else:
        if compareZero:
            logger.TBFail(text, '(abs_prec = %.2e)' % abs_prec)
        else:
            logger.TBFail(text, '(prec = %.2e)' % prec)
