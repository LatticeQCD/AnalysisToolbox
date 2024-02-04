#
# printErrorBars.py
#
# H. Sandmeyer
#
# Methods for printing measurements with error bars as they typically appear in lattice publications. The default is
# 2 significant digits on an error bar in parentheses: measurement = X.XXXX(YY)
#

import numpy as np 
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType


def getValuesFromErrStr(errStr):
    """ Convert A string of the form XX.XX(YY) into a float mean and error bar. Scientific notation not yet supported.

    Args:
        errStr (str): string of the form XX.XXXX(YY).

    Returns:
        float, float: mean, error. 
    """
    checkType(errStr,str)
    try:
        meanStr = errStr.split('(')[0]
        mean    = float(meanStr)
        err     = float(errStr.split('(')[1][:-1])
        if get_exp(err)==0:
            return mean, err
        if not '.' in errStr:
            return mean, err
        dot_index = errStr.find('.')
        parenthesis_index = errStr.find('(')
        N = parenthesis_index - dot_index - 1
        err *= pow(10,-N)
        return mean, err        
    except:
        logger.TBError('Expected string of form XX.XXXX(YY).')


def get_exp(param):
    """ Get the exponent of a number to base 10. """
    result = np.floor( np.log(abs(param))/np.log(10) )
    if np.isnan(result):
        return np.nan
    else:
        return int(result)


def get_err_str(param, param_err, numb_err_dig=2) -> str:
    """ Get the string of a number + error, e.g. 1.234567+-0.324456 --> 12.34(33) (numb_err_dig = 2). """

    checkType(numb_err_dig,int)
    param     = float(param)
    param_err = float(param_err)

    if param_err<=0:
        logger.details('Encountered non-positive error',param_err)
        return param

    if numb_err_dig < 1:
        logger.TBError("Number of error digits has to be larger than 0!")
    if param < 0:
        param = -param
        sign = -1
    else:
        sign = 1

    relnum = get_exp(param_err)
    if np.isnan(relnum):
        return "%.12e" % param

    # index for rounding the error and the parameter
    roundidx = -relnum + numb_err_dig - 1
    paramtmp = param
    param = np.round(param, roundidx)

    # exponent of the actual parameter
    if param == 0.0:
        numdig=get_exp(param_err)
    else:
        numdig = get_exp(param)

    # Number of digits before the dot
    if numdig < 0:
        numdig = 0

    # floor does not support a second index -> we have to multiply and then divide
    param_err *= pow(10, roundidx)
    param_err = np.ceil(param_err) * pow(10, -roundidx)

    # as the exponent might have changed through rounding, we have to recalculate it
    relnum = get_exp(param_err)
    roundidx = -relnum + numb_err_dig - 1
    param = np.round(paramtmp, roundidx)

    # Strings that are shortened later on
    err_str = "%.12lf" % param_err
    param_str = "%.12lf" % param

    # if the error is larger than or equal to 1
    if relnum >= 0:
        # if we need more digits to express the full error than given by numb_err_digits
        if numb_err_dig <= relnum + 1:
            err_str = err_str[:relnum + 1]
            # cut parameter before the dot
            param_str = param_str[:numdig + 1]
        else:
            # if we have to print more digits than relnum
            err_str = err_str[:numb_err_dig + 1]
            # also print after the dot with the relevant number of digits
            param_str = param_str[:numdig + 1 + (numb_err_dig - relnum)]
    else:
        # we don't need the part before the dot here
        err_str = err_str[-relnum + 1:-relnum + numb_err_dig + 1]
        param_str = param_str[:numdig - relnum + numb_err_dig + 1]

    if sign == -1:
        return "-%s(%s)" % (param_str, err_str)
    else:
        return "%s(%s)" % (param_str, err_str)


def get_err_str_exp(param, param_err, exp, numb_err_dig=1, multicon="x") -> str:
    """ Express the number with a exponent. The multiplication icon can be changed. """
    if exp == 0:
        return get_err_str(param, param_err, numb_err_dig)
    param /= pow(10, exp)
    param_err /= pow(10, exp)
    return "%s%s10^%i" % (get_err_str(param, param_err, numb_err_dig), multicon, exp)


def get_err_str_auto(param, param_err, numb_err_dig=1, mulicon="x") -> str:
    """ Automatically express the number with an exponent. The multiplication icon can be changed. """
    exp = get_exp(param)
    exp_err = get_exp(param_err)
    if exp_err > exp:
        exp = exp_err
    if 4 > exp > -4:
        return get_err_str(param, param_err, numb_err_dig)
    return get_err_str_exp(param, param_err, exp, numb_err_dig, mulicon)


def get_err_str_exp_tex(param, param_err, exp, numb_err_dig=1) -> str:
    """ Express the number with an exponent in LaTeX. """
    return "$%s$" % get_err_str_exp(param, param_err, exp, numb_err_dig, "\\cdot ")


def get_err_str_auto_tex(param, param_err, numb_err_dig=1) -> str:
    """ Automatically express the number with an exponent in LaTeX. """
    return "$%s$" % get_err_str_auto(param, param_err, numb_err_dig, "\\cdot ")
