# 
# testing.py                                                               
# 
# D. Clarke 
# 
# Functions to assist with unit testing. 
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import envector
from latqcdtools.math.math import rel_check
from latqcdtools.statistics.statistics import gaudif
from latqcdtools.base.printErrorBars import get_err_str


def print_results(res, res_true, res_err = None, res_err_true = None, text = "", prec = 1e-10, abs_prec = None):
    """ Compares element-by-element the results of res with res_true. (Does the same with res_err and res_err_true,
        if you like.) Carries out with precision prec. Use abs_prec for comparisons with zero. """
    test = True
    compareZero = False

    if abs_prec is not None:
        compareZero = True
    else:
        abs_prec = 1e-14

    res, res_true = envector(res, res_true)
    if res_err is not None:
        res_err = envector(res_err)
    if res_err_true is not None:
        res_err_true = envector(res_err_true)

    for i in range(len(res)):
        if not rel_check(res[i], res_true[i], prec, abs_prec):
            test = False
            logger.info("res[" + str(i) + "] = " + str(res[i]) + " != res_true[" + str(i) + "] = " + str(res_true[i]))
        if res_err is not None and res_err_true is not None:
            if not rel_check(res_err[i], res_err_true[i], prec, abs_prec):
                test = False
                logger.info("res_err[" + str(i) + "] = " + str(res_err[i]) + " != res_err_true[" + str(i) + "] = " + str(res_err_true[i]))

    if test:
        return True
    else:
        if compareZero:
            logger.TBFail(text, '(abs_prec = %.2e)' % abs_prec)
        else:
            logger.TBFail(text, '(prec = %.2e)' % prec)
        return False 


def print_results_iter(res,  res_true, text):
    test = True
    res = np.asarray(res)
    res_true = np.asarray(res_true)
    it = np.nditer(res, flags=['multi_index'])
    while not it.finished:
        if not rel_check(res[it.multi_index], res_true[it.multi_index], 1e-4, 1e-6):
            test = False
            logger.info("res[" + str(it.multi_index) + "] = " + str(res[it.multi_index])
                    + " != res_true[" + str(it.multi_index) + "] = " + str(res_true[it.multi_index]))
        it.iternext()
    if test:
        return True
    else:
        logger.TBFail(text)
        return False


def gaudif_results(res, res_err, res_true, res_err_true, text = "", qcut=0.05):
    """ Compares element-by-element the results of res with res_true using Gaussian difference test, i.e. it checks
        to see whether res and res_true are statistically compatible. """

    test = True

    res, res_true, res_err, res_err_true = envector(res, res_true, res_err, res_err_true)

    for i in range(len(res)):

        q = gaudif(res[i], res_err[i], res_true[i], res_err_true[i])

        if q < qcut:
            test = False
            resstr     = get_err_str(res[i]     ,res_err[i])
            restruestr = get_err_str(res_true[i],res_err_true[i])
            logger.info("res["+str(i)+"] =",resstr,"!= res_true["+str(i)+"] =",restruestr,'[ q =',round(q,2),']')

    if test:
        return True 
    else:
        logger.TBFail(text)
        return False


def concludeTest(lpass):
    if lpass:
        logger.TBPass('All tests passed.')
    else:
        logger.TBError('At least one test failed.')
