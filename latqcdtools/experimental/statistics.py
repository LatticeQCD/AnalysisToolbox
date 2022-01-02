import numpy as np
import mpmath
import latqcdtools.autocorrelation as autocorr
import latqcdtools.tools as tools
import latqcdtools.logger as logger
from latqcdtools.num_deriv import diff_jac










def mean_and_err(data, axis = 0):
    mean = std_mean(data, axis = axis)
    error = std_err(data, axis = axis)
    return mean, error


def mean_and_cov(data, axis = 0):
    mean = std_mean(data, axis = axis)
    cov = calc_cov(data)
    return mean, cov


def mean_and_std_dev(data, axis=0):
    mean = std_mean(data, axis = axis)
    std = std_dev(data, axis = axis)
    return mean, std





def rem_norm_corr(corr, edata):
    """Compute the covariance matrix from a correlation matrix."""
    res = np.zeros((len(corr), len(corr[0])))
    for i in range(len(corr)):
        for j in range(len(corr[0])):
            res[i][j] = corr[i][j] * edata[i]*edata[j]
    return np.array(res)










def aicc_av(aicc, data, errors = None):
    """Weighted average using AICc as weights"""
    data = np.asarray(data)
    errors = np.asarray(errors)
    if errors is None:
        errors = np.ones_like(data)
    aicc = np.asarray(aicc)
    aicc_min = np.min(aicc)
    daicc = aicc - aicc_min
    weights = np.exp(-0.5 * daicc)
    weights /= np.sum(weights)
    return np.sum(weights * data), np.sqrt(np.sum(weights ** 2 * errors ** 2))



def compute_autocorr_impr_err(data, tmax):
    t_arr, corr_arr, corr_err_arr, int_corr_arr, int_corr_err_arr = autocorr.corr(
            data, tmax=tmax, startvalue=0)

    N = len(data)
    err_sqrt = (np.var(data) / N) * (2 * int_corr_arr[tmax - 1])
    return np.sqrt(err_sqrt)


def compute_autocorr_impr_err_arr(data, tmax):
    I=len(data[0])
    K=len(data[0][0])
    err=np.zeros((I, K))
    mean = np.mean(data, axis=0)
    for i in range(len(data[0])):
        print("Column: ",i)
        for k in range(len(data[0][0])):
            err[i,k]=autocorr.compute_autocorr_impr_err(data[:, i, k], tmax)
    return np.transpose(mean), np.transpose(err)

