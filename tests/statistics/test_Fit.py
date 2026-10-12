# 
# test_fit.py                                                               
# 
# H. Sandmeyer, D. Clarke
# 
# Test routines within the fitting module.
# 

import numpy as np
from latqcdtools.statistics.fitting import do_fit, try_fit, Fitter, zipData, unzipData
from latqcdtools.math.math import rel_check
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import timer
from latqcdtools.statistics.statistics import std_mean
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.plotting import plt
from latqcdtools.base.readWrite import readTable


SHOWPLOTS = False


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
        res += np.array([(data[:, l] - mean)]).transpose().dot(np.array([(data[:, l] - mean)]))
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

    # First we test zip and unzip for 2d and 3d. Hopefully this gives the user a sense
    # of how the indexing works.

    xdata = np.array([1,2,3,4])
    ydata = np.array([5,6,7])
    zdata = np.array([8,9])

    xtest, ytest = unzipData(zipData(xdata,ydata))

    xREF = np.array([1,1,1,2,2,2,3,3,3,4,4,4])
    yREF = np.array([5,6,7,5,6,7,5,6,7,5,6,7])

    lpass *= print_results(xtest,xREF,text='zip-unzip x 2d')
    lpass *= print_results(ytest,yREF,text='zip-unzip y 2d')

    xtest, ytest, ztest = unzipData(zipData(xdata,ydata,zdata))

    xREF = np.array([1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4])
    yREF = np.array([5,5,6,6,7,7,5,5,6,6,7,7,5,5,6,6,7,7,5,5,6,6,7,7])
    zREF = np.array([8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9,8,9])

    lpass *= print_results(xtest,xREF,text='zip-unzip x 3d')
    lpass *= print_results(ytest,yREF,text='zip-unzip y 3d')
    lpass *= print_results(ztest,zREF,text='zip-unzip z 3d')



    logger.info("Testing quadratic fit with expansion of parameters...")

    xdata, ydata, edata = readTable("../../datasets/wurf.dat", usecols=(0,2,3))

    res_true     = [-1.930355e+00, 6.747380e+00, -6.979050e-02]
    res_err_true = [ 9.072495e-02, 3.357190e-01,  2.676424e-01]

    start_params=np.array([1,1,1])

    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, start_params, algorithm="TNC", grad = grad_fit_func, norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact TNC",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, start_params, algorithm="TNC", grad = grad_fit_func, 
                             norm_err_chi2=True,error_strat='hessian')
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact TNC, hessian error strat",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, start_params, algorithm="L-BFGS-B", derive_chisq = True, norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true,"Numerical L-BFGS-B using built-in derivative",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, start_params, algorithm="SLSQP", derive_chisq= True, norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true,"Numerical SLSQP using built-in derivative",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, start_params, algorithm="Powell", norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true, "Powell quadratic ",prec=EPSILON)


    res, res_err, _ = do_fit(fit_func, xdata, ydata, edata, start_params, algorithm="nonlin", norm_err_chi2=True) 
    lpass *= print_results(res, res_true, res_err, res_err_true, "Nonlinear least-squares ",prec=EPSILON)



    logger.info("Testing correlator fit...")

    xdata, ydata, edata = readTable("../../datasets/corr.dat", usecols=(0,1,2))

    res_true = [5.088129e-05, 2.943403e-01]
    res_err_true = [5.042611e-08, 8.380914e-05]

    start_params=np.array([1,1])

    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, start_params, grad=grad_one_state,
                             args=(64,), norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact curve_fit",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, np.diag(edata)**2, start_params, grad=grad_one_state,
                             args=(64,), norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true, "Diagonal correlation matrix",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, start_params, args=(64,), use_diff = False,
                             norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true, "Numerical curve_fit with difference quotient applied on chisquare",
                           prec=EPSILON)


    # Numerical derivative gives a slightly different result
    res_err_true = [5.0425819803e-08, 8.38114689761e-05]
    res, res_err, _ = do_fit(one_state, xdata, ydata, edata, start_params, args=(64,), use_diff = True,
                             norm_err_chi2=True, algorithm="curve_fit")
    lpass *= print_results(res, res_true, res_err, res_err_true,"Numerical curve_fit with difference quotient",prec=EPSILON)


    xdata, data, nconfs = readCorrelatorTable("../../datasets/corr_pure.dat",1,2)

    cov   = calc_cov_OLD(data)
    ydata = std_mean(data, axis = 1)
    cov_true = readTable("../../datasets/cov.txt")
    cov_test = True
    for i in range(len(cov)):
        for j in range(len(cov[i])):
            if not rel_check(cov[i, j], cov_true[i, j]):
                cov_test = False
                logger.info(f"cov[{i},{j}] = {cov[i, j]} != cov_true[{i},{j}] = {cov_true[i, j]}")
    if not cov_test:
        logger.TBFail("Covariance matrix test")
        lpass = False

    res, res_err, _ = do_fit(one_state, xdata, ydata, cov / nconfs, np.array(res_true), grad=grad_one_state,
                             args=(64,), norm_err_chi2=True,algorithm="curve_fit")
    res_true = [4.988713e-05, 2.950030e-01]
    res_err_true = [1.176005e-06, 5.573209e-04]
    lpass *= print_results(res, res_true, res_err, res_err_true, "Exact curve_fit for correlated data",prec=EPSILON)


    res, res_err, _ = do_fit(one_state, xdata, ydata, cov / nconfs, np.array(res_true), args=(64,),
                             algorithm = "curve_fit", norm_err_chi2=True)
    lpass *= print_results(res, res_true, res_err, res_err_true, "Numerical curve_fit for correlated data",prec=EPSILON)



    logger.info("Testing Bayesian fit...")

    prior        = np.array([5e-05,3e-01])
    priorsigma   = np.array([5e-05,3e-01])
    # lsqfit 13.2.3 errors times sqrt(chi^2/dof), since norm_err_chi2=True.
    res_err_true = [1.148828061442774e-06, 0.0005444934695280138]
    res, res_err, _, stats = try_fit(one_state, xdata, ydata, cov / nconfs, priorval = prior, priorsigma = priorsigma, 
                                     args=(64,), norm_err_chi2=True, detailedInfo=True,
                                     algorithms = ["L-BFGS-B", "TNC", "Powell" ,"Nelder-Mead", "dogleg", "trust-ncg"])
    lpass *= print_results(res, res_true, res_err, res_err_true, "Constraint fit",prec=EPSILON)
    res, res_err, _ = try_fit(one_state, xdata, ydata, cov / nconfs, priorval = prior, priorsigma = priorsigma, 
                              args=(64,), norm_err_chi2=True, algorithms = ["nonlin"])
    lpass *= print_results(res, res_true, res_err, res_err_true, "Constraint fit, nonlin",prec=EPSILON)

    # lsqfit 13.2.3 gives 542.4454270481813 for this fit.
    lpass *= print_results(stats['logGBF'], 542.4454270418508     , text='logGBF', prec=1e-9)
    lpass *= print_results(stats['chi2']  , 95.21248701567058     , text='chi2')
    lpass *= print_results(stats['BAIC']  , 99.21220447790006     , text='BAIC')
    lpass *= print_results(stats['AIC']   , -1102.839241438541    , text='AIC')
    lpass *= print_results(stats['AICc']  , -1102.5465585117117   , text='AICc')
    lpass *= print_results(stats['Q']     , 1.2269138243695367e-05, text='Q')

    # A fit with more parameters than data, so the priors determine the errors of the higher coefficients.
    # Reference values are from lsqfit 13.2.3.
    def quartic(x, p):
        return p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3 + p[4]*x**4
    res_true     = [0.96554345, 1.38433448, 0.15002567, 0.04753187, 0.0143223]
    res_err_true = [0.07654904, 0.47778756, 0.92638066, 0.98949346, 0.99866696]
    for strat in ['propagation','hessian']:
        res, res_err, _ = try_fit(quartic, np.array([0.1,0.2,0.3]), np.array([1.1,1.25,1.4]), np.array([0.05,0.05,0.05]),
                                  priorval=np.array([1.,1.,0.,0.,0.]), priorsigma=np.ones(5), error_strat=strat,
                                  algorithms=["L-BFGS-B", "TNC", "Powell" ,"Nelder-Mead"])
        lpass *= print_results(res, res_true, res_err, res_err_true, "Prior-dominated fit, "+strat, prec=1e-6)
    res, res_err, _ = try_fit(quartic, np.array([0.1,0.2,0.3]), np.array([1.1,1.25,1.4]), np.array([0.05,0.05,0.05]),
                              priorval=np.array([1.,1.,0.,0.,0.]), priorsigma=np.ones(5), algorithms=["nonlin"])
    lpass *= print_results(res, res_true, res_err, res_err_true, "Prior-dominated fit, nonlin", prec=1e-6)

    # Uncorrelated line fit whose errors span seven orders of magnitude. The SVD cut acts on the correlation
    # matrix, so it must leave these errors alone, and the fit errors must equal the exact weighted least squares.
    xline = np.array([0.,1.,2.,3.])
    yline = np.array([1.1,2.9,5.0,7.2])
    eline = np.array([1e-7,1.,1.,1.])
    J     = np.array([np.ones(4), xline]).T
    pcov  = np.linalg.inv(J.T @ np.diag(1/eline**2) @ J)
    def line(x, p):
        return p[0] + p[1]*x
    res, res_err, _ = try_fit(line, xline, yline, eline, start_params=[0,0], algorithms=["curve_fit"])
    lpass *= print_results(res_err, np.sqrt(np.diag(pcov)), text="Errors spanning many orders of magnitude", prec=1e-8)

    # Only p[0]+p[1] is constrained, so pcov does not exist. The fit must fail rather than report finite errors.
    def degenerate(x, p):
        return (p[0]+p[1])*x
    def grad_degenerate(x, p):
        return np.array([x, x])
    try:
        try_fit(degenerate, xline, yline, np.ones(4), grad=grad_degenerate, start_params=[1,1], algorithms=["curve_fit"])
        logger.TBFail('Degenerate fit did not fail')
        lpass = False
    except Exception:
        logger.TBPass('Degenerate fit fails')

    timey.printTiming()


    logger.info("Testing fit plots...")

    xdata, ydata, edata = readTable("../../datasets/wurf.dat", usecols=(0,2,3))

    fitter = Fitter(fit_func,xdata,ydata,edata, norm_err_chi2=True)
    fitter.do_fit(start_params=np.array([1,1,1]))

    fitter.plot_fit()
    if SHOWPLOTS:
        plt.show()


    concludeTest(lpass)


if __name__ == '__main__':
    testFit()
