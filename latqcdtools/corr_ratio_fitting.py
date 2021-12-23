from latqcdtools.fitting import print_res, print_scl, Fitter
from latqcdtools.corr_fitting import CorrFitter
from latqcdtools.statistics import mean_and_err, mean_and_std_dev
import latqcdtools.logger as logger
import numpy as np



def corr_ratio_err(data):
    mean, err=mean_and_err(data, axis = 1)
    res = [ mean[i] / mean[i+1] for i in range(len(data) - 1 )]
    res_err = [ np.sqrt( (err[i] / mean[i + 1])**2 +
        (err[i+1] / mean[i + 1]**2)**2) for i in range(len(data) - 1) ]
    return res, res_err

def corr_ratio_std(data):
    mean, std = mean_and_std_dev(data, axis = 1)
    res = [ mean[i] / mean[i+1] for i in range(len(data) - 1 )]
    res_std = [ np.sqrt( (std[i] / mean[i + 1])**2
        + (std[i+1] * mean[i] / mean[i + 1]**2)**2) for i in range(len(data) - 1) ]
    return res, res_std

def corr_ratio(data):
    mean=np.mean(data, axis = 1)
    return [ mean[i] / mean[i+1] for i in range(len(mean) -1 )]

def corr_ratio_direct(mean):
    return [ mean[i] / mean[i+1] for i in range(len(mean) -1 )]



class CorrFitterRatio:

    def __init__(self, xdata, ydata, edata, xdata_est, ydata_est, edata_est, nstates = 1, Nt = None):


        if Nt is None:
            self._Nt = len(ydata) + 1
        else:
            self._Nt = Nt

        self._estimator = CorrFitter(xdata_est, ydata_est, edata_est,
                nstates = nstates, nstates_osc = 0, Nt = self._Nt)

        if nstates == 1:
            self._fitter = Fitter(one_state_ratio, xdata, ydata, edata,
                    grad = grad_one_state_ratio, hess = hess_one_state_ratio, 
                    args = (self._Nt,))

        elif nstates == 2:
            self._fitter = Fitter(two_state_ratio, xdata, ydata, edata,
                    grad = grad_two_state_ratio, hess = hess_two_state_ratio, 
                    args = (self._Nt,))
        else:
            raise ValueError("Maximal two states for ratio fits")
        
        self._nstates = nstates
        self._xdata = np.asarray(xdata)


    def plot_eig(self, **kwargs):
        return self._fitter.plot_eig(**kwargs)

    def plot_cov(self, **kwargs):
        return self._fitter.plot_cov(**kwargs)


    def _swap(self, res, res_err):
        if (abs(res[1]) > abs(res[2]) > 1e-5
                or (abs(res[1]) < 1e-5 < abs(res[2]))):

            res[1], res[2] = abs(res[2]), abs(res[1]),
            res_err[1], res_err[2] = abs(res_err[2]), abs(res_err[1]),



    def corr_fit(self, xmin = -np.inf, xmax = np.inf, start_params = None,
            priorsigma = None, priorval = None, nstates = 1):
        if priorval is not None:
            start_params = np.copy(priorval)
            priorval = np.copy(priorval)
            priorsigma = np.copy(priorsigma)

        if start_params is None:
            try:
                tmp_params = self._estimator.corr_fit(xmin,
                        xmax, correlated = False)[0]
                if self._nstates == 2:
                    start_params = [tmp_params[0]/tmp_params[2], tmp_params[1], tmp_params[3]]
                else:
                    start_params = [tmp_params[1]]
            except Exception as e:
                logger.warn(e)
                logger.warn("\nEstimating start parameters for ratio fit failed. Try direct fit\n")

        print_res("Final starting parameters for ratio fit", start_params, level = "INFO")
        res, res_err, chi_dof, pcov = self._fitter.try_fit(xmin = xmin, xmax = xmax,
                start_params = start_params, priorval = priorval, priorsigma = priorsigma,
                ret_pcov = True)

        aicc = self._fitter.aicc(res, xmin, xmax, correlated = None)
        if self._nstates > 1:
            self._swap(res, res_err)

        print_res("Fit result from ratio fit", res, res_err, chi_dof, level = "INFO")
        print_scl("AICc", aicc, level = "INFO")
        return res, res_err, chi_dof, aicc, pcov


    def plot_corr(self, filename, params, params_err, ranges, xlabel = None, ylabel = None,
                 title=None, notex=False, plot_ylog=False, xmin = 1, **kwargs):


        if xlabel is None:
            xlabel = "$n_{\\tau/\\sigma}$"

        if ylabel is None:
            ylabel = "$G(n_{\\tau/\\sigma})/G(n_{\\tau/\\sigma} + 1)$"

        if np.max(ranges) <= self._Nt/2:
            xmax = self._Nt/2
        else:
            xmax = np.max(self._xdata)


        self._fitter.plot_fit(filename, params, params_err, notex = notex, ylog = plot_ylog,
                args_data = {'xmin': xmin, 'xmax' : xmax, 'xlabel' : xlabel, 'ylabel' : ylabel},
                ranges = ranges, fix_ylim = True, **kwargs)




def one_state_ratio(x, m, Nt):
    return np.cosh(m * (Nt / 2 - x)) / np.cosh(m * (Nt/2 - x - 1))

def two_state_ratio(x, A2, m2, m1, Nt):
    return (np.cosh(m1 * (Nt / 2 - x)) + A2 * np.cosh(m2 * (Nt / 2 - x)))\
            / (np.cosh(m1 * (Nt / 2 - x - 1)) + A2 * np.cosh(m2 * (Nt / 2 - x - 1)))

def grad_one_state_ratio(x, m, Nt):
    if x == Nt / 2 - 1:
        return np.array([np.sinh(m * (Nt / 2 - x)) * (Nt / 2 - x) ])
    return np.array([np.sinh(m * (Nt / 2 - x)) * (Nt / 2 - x)
        / np.cosh(m * (Nt/2 - x - 1)) + np.cosh(m * (Nt / 2 - x))
        / np.sinh(m * (Nt/2 - x - 1))**2 * (Nt/2 - x - 1)])

def hess_one_state_ratio(x, m, Nt):
    if x == Nt / 2 - 1:
        return np.array([[np.cosh(m * (Nt / 2 - x)) * (Nt / 2 - x)**2 ]])
    return np.array([[0.25 * 1/np.sinh(m * (1 - Nt/2 + x)) *
        ((Nt - 2 * x)**2 * np.cosh(0.5 * m * (Nt - 2 * x))
        - (Nt - 2 * (1 + x))**2 * np.cosh(0.5 * m * (Nt - 2 * x))
        * 1/np.sinh(m * (1 - Nt/2 + x))**2 - 2 * (Nt - 2 * x) * (Nt - 2 * (1 + x)) * 
        np.sinh(0.5 * m * (Nt - 2 * x)) * np.tanh(0.5 * m * (Nt - 2 * (1 + x)))
        + (Nt - 2 * (1 + x))**2 * np.cosh(0.5 * m * (Nt - 2 * x))
        * np.tanh(0.5 * m * (Nt - 2 * (1 + x)))**2)]])


def grad_two_state_ratio(x, A2, m2, m1, Nt):
    if x == Nt / 2 - 1:
        return np.array([
            np.cosh(m2 * (Nt / 2 - x)),
            A2 * np.sinh(m2 * (Nt / 2 - x)) * (Nt / 2 - x),
            np.sinh(m1 * (Nt / 2 - x)) * (Nt / 2 - x)])
    m1c=np.cosh(m1 * (Nt/2 - x))
    m1s=np.sinh(m1 * (Nt/2 - x))
    m2c=np.cosh(m2 * (Nt/2 - x))
    m2s=np.sinh(m2 * (Nt/2 - x))
    m1c_1=np.cosh(m1 * (-1 + Nt/2 - x))
    m1s_1=np.sinh(m1 * (-1 + Nt/2 - x))
    m2c_1=np.cosh(m2 * (-1 + Nt/2 - x))
    m2s_1=np.sinh(m2 * (-1 + Nt/2 - x))
    ntx = Nt/2 - x
    ntx1 = Nt/2 - x - 1

    res = np.array([m2c/(m1c_1 + A2 *m2c_1) -
        (m2c_1 * (m1c + A2 * m2c)) /(m1c_1 + A2 * m2c_1)**2,
         -((A2 * ntx1 * (m1c + A2 * m2c) * m2s_1) /(m1c_1 + A2 * m2c_1)**2)
         + (A2 * ntx * m2s)/(m1c_1 + A2 * m2c_1),
        -((ntx1 * (m1c + A2 * m2c) * m1s_1) / (m1c_1 + A2 * m2c_1)**2)
        + (ntx * np.sinh(m1 * ntx))/(m1c_1 + A2 * m2c_1)])
    return res




# Hopefully mathematica is correct. What a mess...
def hess_two_state_ratio(x, A2, m2, m1, Nt):
    if x == Nt / 2 - 1:
        return np.array([[0, np.sinh(m2 * (Nt / 2 - x)) * (Nt / 2 - x), 0],
                [np.sinh(m2 * (Nt / 2 - x)) * (Nt / 2 - x), A2
                    * np.cosh(m2 * (Nt / 2 - x)) * (Nt / 2 - x)**2, 0],
                [0, 0, np.cosh(m1 * (Nt / 2 - x)) * (Nt / 2 - x)**2]])

    m1c=np.cosh(m1 * (Nt/2 - x))
    m1s=np.sinh(m1 * (Nt/2 - x))
    m2c=np.cosh(m2 * (Nt/2 - x))
    m2s=np.sinh(m2 * (Nt/2 - x))
    m1c_1=np.cosh(m1 * (-1 + Nt/2 - x))
    m1s_1=np.sinh(m1 * (-1 + Nt/2 - x))
    m2c_1=np.cosh(m2 * (-1 + Nt/2 - x))
    m2s_1=np.sinh(m2 * (-1 + Nt/2 - x))
    ntx = Nt/2 - x
    ntx1 = Nt/2 - x - 1
    res = np.array([
            [
            -(( 2 * m2c_1 * m2c)/(m1c_1 +  A2*m2c_1)**2)
            + ( 2 * m2c_1**2 * (m1c +  A2*m2c))/(m1c_1 +  A2*m2c_1)**3,
            -(( A2 * ntx1 * m2c * m2s_1)/(m1c_1 +  A2*m2c_1)**2)
            + (2 * A2 * ntx1 * m2c_1 * (m1c +  A2*m2c) * m2s_1)/(m1c_1 + A2*m2c_1)**3
            - (ntx1 * (m1c + A2*m2c) * m2s_1)/(m1c_1 +  A2*m2c_1)**2
            - ( A2 * ntx * m2c_1 * m2s)/(m1c_1 +  A2*m2c_1)**2 + (ntx * m2s) \
                    /( m1c_1 +  A2*m2c_1),
            -((ntx1 * m2c * m1s_1)/(m1c_1 +  A2*m2c_1)**2)
            + (2 * ntx1 * m2c_1 * (m1c +  A2*m2c) * m1s_1)/(m1c_1 + A2*m2c_1)**3
            - (ntx * m2c_1 * m1s)/(m1c_1 +  A2*m2c_1)**2
            ],[
            -(( A2 * ntx1 * m2c * m2s_1)/(m1c_1 +  A2*m2c_1)**2)
            + (2 * A2 * ntx1 * m2c_1 * (m1c +  A2*m2c) * m2s_1)/(m1c_1 + A2*m2c_1)**3
            - (ntx1 * (m1c + A2*m2c) * m2s_1)/(m1c_1 +  A2*m2c_1)**2
            - ( A2 * ntx * m2c_1 * m2s)/(m1c_1 +  A2*m2c_1)**2 + (ntx * m2s)\
                    / ( m1c_1 + A2*m2c_1),
            ( A2 * ntx**2 * m2c)/( m1c_1 + A2*m2c_1)
            - ( A2 * ntx1**2 * m2c_1 * (m1c +  A2*m2c))/(m1c_1 +  A2*m2c_1)**2
            + ( 2 * A2**2 * ntx1**2 * (m1c +  A2*m2c) * m2s_1**2)/(m1c_1 +  A2*m2c_1)**3
            - ( 2 * A2**2 * ntx1 * ntx * m2s_1 * m2s)/(m1c_1 +  A2*m2c_1)**2,
            (2 * A2 * ntx1**2 * (m1c + A2*m2c) * m1s_1 * m2s_1)/(m1c_1 + A2*m2c_1)**3
            - ( A2 * ntx1 * ntx * m2s_1 * m1s)/(m1c_1 +  A2*m2c_1)**2
            - ( A2 * ntx1 * ntx * m1s_1 * m2s)/(m1c_1 +  A2*m2c_1)**2
            ],[
            -((ntx1 * m2c * m1s_1)/(m1c_1 +  A2*m2c_1)**2)
            + (2 * ntx1 * m2c_1 * (m1c +  A2*m2c) * m1s_1)/(m1c_1 + A2*m2c_1)**3
            - (ntx * m2c_1 * m1s)/(m1c_1 +  A2*m2c_1)**2,
            (2 * A2 * ntx1**2 * (m1c + A2*m2c) * m1s_1 * m2s_1)/(m1c_1 + A2*m2c_1)**3
            - ( A2 * ntx1 * ntx * m2s_1 * m1s)/(m1c_1 +  A2*m2c_1)**2
            - ( A2 * ntx1 * ntx * m1s_1 * m2s)/(m1c_1 +  A2*m2c_1)**2,
                (ntx**2 * m1c)/( m1c_1 + A2*m2c_1)
                - (ntx1**2 * m1c_1 * (m1c + A2*m2c))/(m1c_1 + A2*m2c_1)**2
            + ( 2 * ntx1**2 * (m1c +  A2*m2c) * m1s_1**2)/(m1c_1 +  A2*m2c_1)**3
            - ( 2 * ntx1 * ntx * m1s_1 * m1s)/(m1c_1 +  A2*m2c_1)**2
            ]] ) 

    return res