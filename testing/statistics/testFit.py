# 
# test_fit.py                                                               
# 
# H. Sandmeyer, D. Clarke
# 
# Test routines within the fitting module.
# 

import numpy as np
from latqcdtools.statistics.fitting import do_fit, try_fit, Fitter
from latqcdtools.math.math import print_results, rel_check
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import timer
from latqcdtools.statistics.statistics import std_mean
from latqcdtools.base.readWrite import readTable, readCorrelatorTable
from latqcdtools.base.plotting import clearPlot, latexify, plt
from latqcdtools.base.readWrite import readTable


logger.set_log_level('INFO')
TESTPLOTS = False


''' Quadratic fit '''
def fit_func(x, a, b, c):
    return a * x**2 + b * x + c

def grad_fit_func(x, a, b, c):
    return np.array([x**2, x, np.ones(len(x))])

def hess_fit_func(x, a, b, c):
    return np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])


''' Correlator fit. '''
def one_state(x, A, m, Nt):
    return A * np.cosh(m * (Nt / 2 - x))

def grad_one_state(x, A, m, Nt):
    return np.array([np.cosh(m * (Nt / 2 - x)), A * np.sinh(m * (Nt / 2 - x)) * (Nt / 2 - x)])

def hess_one_state(x, A, m, Nt):
    return np.array([ [0, np.sinh(m * (Nt / 2 - x)) * (Nt / 2 - x)],
                      [np.sinh(m*(Nt/2-x))*(Nt/2-x), A*np.cosh(m*(Nt/2-x)) * (Nt/2-x)**2] ])

def calc_cov_OLD(data):
    data = np.asarray(data)
    mean = np.mean(data, axis=1)
    res = np.zeros((len(data), len(data)))
    for l in range(0, len(data[0])):
        res += np.array([(data[:, l] - mean)]).transpose().dot(
            np.array([(data[:, l] - mean)]))
    return 1 / (len(data[0]) - 1) * res


EPSILON=1e-4


def testFit():

    timey = timer()

    logger.info("Testing quadradic fit with expansion of parameters...")

    xdata, ydata, edata = readTable("wurf.dat", usecols=(0,2,3))

    res_true     = [-1.930355e+00, 6.747380e+00, -6.979050e-02]
    res_err_true = [ 9.072495e-02, 3.357190e-01,  2.676424e-01]


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="TNC", grad = grad_fit_func, norm_err_chi2=True)
    print_results(res, res_true, res_err, res_err_true, "Exact TNC",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="L-BFGS-B", derive_chisq = True, norm_err_chi2=True)
    print_results(res, res_true, res_err, res_err_true,"Numerical L-BFGS-B using built-in derivative",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="SLSQP", derive_chisq= True, norm_err_chi2=True)
    print_results(res, res_true, res_err, res_err_true,"Numerical SLSQP using built-in derivative",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, [1, 1, 1], algorithm="Powell", norm_err_chi2=True)
    print_results(res, res_true, res_err, res_err_true, "Powell quadratic ",prec=EPSILON)



    logger.info("Testing correlator fit...")


    xdata, ydata, edata = readTable("corr.dat", usecols=(0,1,2))


    res_true = [5.088129e-05, 2.943403e-01]
    res_err_true = [5.042611e-08, 8.380914e-05]


    res, res_err, _ = do_fit(one_state, xdata, ydata, edata,[1, 1], grad=grad_one_state, hess=hess_one_state,
                                     args=(64,), norm_err_chi2=True, algorithm="curve_fit")
    print_results(res, res_true, res_err, res_err_true, "Exact curve_fit",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, np.diag(edata)**2, [1, 1], grad=grad_one_state,
                                     hess=hess_one_state, args=(64,), norm_err_chi2=True,
                                     algorithm="curve_fit")
    print_results(res, res_true, res_err, res_err_true, "Diagonal correlation matrix",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, [1, 1], args=(64,), use_diff = False,
                                     norm_err_chi2=True, algorithm="curve_fit")
    print_results(res, res_true, res_err, res_err_true, "Numerical curve_fit with difference quotient applied on chisquare"
                  ,prec=EPSILON)


    # Numerical derivative gives a slightly different result
    res_err_true = [5.0425819803e-08, 8.38114689761e-05]
    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, [1, 1], args=(64,), use_diff = True,
                                     norm_err_chi2=True, algorithm="curve_fit", )
    print_results(res, res_true, res_err, res_err_true,"Numerical curve_fit with difference quotient",prec=EPSILON)


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
    if cov_test:
        logger.TBPass("Covariance matrix test")
    else:
        logger.TBFail("Covariance matrix test")


    res, res_err, _ = do_fit(one_state, xdata, ydata, cov / nconfs, res_true, grad=grad_one_state,
                                     hess=hess_one_state, args=(64,), norm_err_chi2=True,
                                     algorithm="curve_fit")
    res_true = [4.988713e-05, 2.950030e-01]
    res_err_true = [1.176005e-06, 5.573209e-04]
    print_results(res, res_true, res_err, res_err_true, "Exact curve_fit for correlated data",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, cov / nconfs, res_true, args=(64,),
                                     algorithm = "curve_fit", norm_err_chi2=True)
    print_results(res, res_true, res_err, res_err_true, "Numerical curve_fit for correlated data",prec=EPSILON)



    logger.info("Testing Bayesian fit...")

    prior        = [5e-05,3e-01]
    prior_err    = [5e-05,3e-01]
    res_err_true = [1.148967933055725e-06, 0.0005445077355994696]
    res, res_err, _, _, _ = try_fit(one_state, xdata, ydata, cov / nconfs, priorval = prior, priorsigma = prior_err, args=(64,),
                                    norm_err_chi2=True, detailedInfo=True,
                                    algorithms = ["L-BFGS-B", "TNC", "Powell" ,"Nelder-Mead", "dogleg", "trust-ncg"])
    print_results(res, res_true, res_err, res_err_true, "Constraint fit",prec=EPSILON)

    timey.printTiming()


    if TESTPLOTS:

        logger.info("Testing fit plots...")
        latexify()

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

    logger.TBPass("No problems encountered.")

if __name__ == '__main__':
    testFit()
