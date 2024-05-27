# 
# test_fit.py                                                               
# 
# H. Sandmeyer, D. Clarke
# 
# Test routines within the fitting module.
# 

import numpy as np
from latqcdtools.statistics.fitting import do_fit, try_fit, Fitter
from latqcdtools.math.math import rel_check
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import timer
from latqcdtools.statistics.statistics import std_mean
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.plotting import clearPlot, plt
from latqcdtools.base.readWrite import readTable


logger.set_log_level('INFO')
TESTPLOTS = False


''' Quadratic fit '''
def fit_func(x, p):
    a, b, c = p
    return a * x**2 + b * x + c

def grad_fit_func(x, p):
    return np.array([x**2, x, np.ones(len(x))])


''' Correlator fit. '''
def one_state(x, p, Nt):
    A, m = p
    return A * np.cosh(m * (Nt / 2 - x))

def grad_one_state(x, p, Nt):
    A, m = p
    return np.array([np.cosh(m * (Nt / 2 - x)), A * np.sinh(m * (Nt / 2 - x)) * (Nt / 2 - x)])


def calc_cov_OLD(data):
    data = np.asarray(data)
    mean = np.mean(data, axis=1)
    res = np.zeros((len(data), len(data)))
    for l in range(0, len(data[0])):
        res += np.array([(data[:, l] - mean)]).transpose().dot(
            np.array([(data[:, l] - mean)]))
    return 1 / (len(data[0]) - 1) * res


def readCorrelatorTable(filename, col1=1, col2=2, symmetrize = False):
    """ Read a table organized like so:
    t=0    meas    meas_err conf
    t=1    meas    meas_err conf
    ...
    t=Nt/2 meas    meas_err conf
    t=0    meas    meas_err conf+1
    ...
    """
    try:
        # To support input file streams
        ins = open(filename, "r")
        close = True
    except TypeError:
        ins = filename
        close = False
    data_dict = {}
    for line in ins:
        if line.startswith('#') or len(line) < 2:
            continue
        lineelems = line.strip().split()
        try:
            Nt = int(lineelems[col1 - 1])
        except ValueError:
            Nt = float(lineelems[col1 - 1])
        corr = float(lineelems[col2 - 1])
        if Nt not in data_dict:
            data_dict[Nt] = []
        data_dict[Nt].append(corr)
    xdata = list(sorted(data_dict))
    data = [ data_dict[key] for key in sorted(data_dict.keys()) ]
    Nt = len(data)
    if symmetrize:
        if max(xdata) != Nt - 1:
            logger.TBError("The number of x values does not correspond to the largest of its values")
        if Nt % 2 != 0:
            logger.TBError("Nt must be even!")

        for i in range(len(data[0])):
            for nt in range(1, int(len(data)/2)):
                data[nt][i] = data[Nt - nt][i] = (data[nt][i] + data[Nt - nt][i]) / 2
    if close:
        ins.close()
    return np.array(xdata), np.array(data), len(data[0])


EPSILON=1e-4


def testFit():

    timey = timer()
    lpass = True

    logger.info("Testing quadradic fit with expansion of parameters...")

    xdata, ydata, edata = readTable("wurf.dat", usecols=(0,2,3))

    res_true     = [-1.930355e+00, 6.747380e+00, -6.979050e-02]
    res_err_true = [ 9.072495e-02, 3.357190e-01,  2.676424e-01]


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="TNC", grad = grad_fit_func, norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact TNC",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="TNC", grad = grad_fit_func, 
                             norm_err_chi2=True,error_strat='hessian')
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact TNC, hessian error strat",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="L-BFGS-B", derive_chisq = True, norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true,"Numerical L-BFGS-B using built-in derivative",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="SLSQP", derive_chisq= True, norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true,"Numerical SLSQP using built-in derivative",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="Powell", norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true, "Powell quadratic ",prec=EPSILON)



    logger.info("Testing correlator fit...")


    xdata, ydata, edata = readTable("corr.dat", usecols=(0,1,2))


    res_true = [5.088129e-05, 2.943403e-01]
    res_err_true = [5.042611e-08, 8.380914e-05]


    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, [1, 1], grad=grad_one_state,
                             args=(64,), norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact curve_fit",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, np.diag(edata)**2, [1, 1], grad=grad_one_state,
                             args=(64,), norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true, "Diagonal correlation matrix",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, [1, 1], args=(64,), use_diff = False,
                             norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true, "Numerical curve_fit with difference quotient applied on chisquare",
                           prec=EPSILON)


    # Numerical derivative gives a slightly different result
    res_err_true = [5.0425819803e-08, 8.38114689761e-05]
    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, [1, 1], args=(64,), use_diff = True,
                             norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true,"Numerical curve_fit with difference quotient",prec=EPSILON)


    xdata, data, nconfs = readCorrelatorTable("corr_pure.dat",1,2)

    cov   = calc_cov_OLD(data)
    ydata = std_mean(data, axis = 1)
    cov_true = readTable("cov.txt")
    cov_test = True
    for i in range(len(cov)):
        for j in range(len(cov[i])):
            if not rel_check(cov[i, j], cov_true[i, j]):
                cov_test = False
                logger.info("cov[" + str(i) + "," + str(j) + "] = " + str(cov[i, j])
                        + " != cov_true[" + str(i) + "," + str(j) + "] = " + str(cov_true[i, j]))
    if not cov_test:
        logger.TBFail("Covariance matrix test")
        lpass = False

    res, res_err, _ = do_fit(one_state, xdata, ydata, cov / nconfs, res_true, grad=grad_one_state,
                             args=(64,), norm_err_chi2=True,algorithm="curve_fit")
    res_true = [4.988713e-05, 2.950030e-01]
    res_err_true = [1.176005e-06, 5.573209e-04]
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact curve_fit for correlated data",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, cov / nconfs, res_true, args=(64,),
                             algorithm = "curve_fit", norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true, "Numerical curve_fit for correlated data",prec=EPSILON)



    logger.info("Testing Bayesian fit...")

    prior        = [5e-05,3e-01]
    priorsigma   = [5e-05,3e-01]
    res_err_true = [1.148967933055725e-06, 0.0005445077355994696]
    res, res_err, _, stats = try_fit(one_state, xdata, ydata, cov / nconfs, priorval = prior, priorsigma = priorsigma, 
                                     args=(64,), norm_err_chi2=True, detailedInfo=True,
                                     algorithms = ["L-BFGS-B", "TNC", "Powell" ,"Nelder-Mead", "dogleg", "trust-ncg"])
    lpass *= print_results(res, res_true, res_err, res_err_true, "Constraint fit",prec=EPSILON)

    lpass *= print_results(stats['logGBF'], 542.3120190935467     , text='logGBF')
    lpass *= print_results(stats['chi2']  , 95.21248701567058     , text='chi2')
    lpass *= print_results(stats['BAIC']  , 99.21220447790006     , text='BAIC')
    lpass *= print_results(stats['AIC']   , -1080.6240381870934   , text='AIC')
    lpass *= print_results(stats['AICc']  , -1080.3573715204268   , text='AICc')
    lpass *= print_results(stats['Q']     , 1.2269138243695367e-05, text='Q')

    timey.printTiming()


    if TESTPLOTS:

        logger.info("Testing fit plots...")

        xdata, ydata, edata = readTable("wurf.dat", usecols=(0,2,3))

        fitter = Fitter(fit_func,xdata,ydata,edata, norm_err_chi2=True)
        fitter.do_fit(start_params=None)

        fitter.plot_fit()
        plt.show()
        clearPlot()
        fitter.plot_cor()
        plt.show()
        clearPlot()
        fitter.plot_eig()
        plt.show()


    concludeTest(lpass)


if __name__ == '__main__':
    testFit()
