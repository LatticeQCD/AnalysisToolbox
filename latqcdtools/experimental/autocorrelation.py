import numpy as np


def corr(data, tmax=100, startvalue=0):
    """Error computation according to Box, Jenkins, Reinsel: "Time series analysis", third edition, p. 33, equation
    2.1.13. This method of calculating autocorrelation does not remove the bias."""

    int_corr = 0.5
    sq_sum = 0.0

    t_arr = []
    corr_arr = []
    int_corr_arr = []
    corr_err_arr = []
    int_corr_err_arr = []

    mean = np.mean(data,axis=0)
    norm = np.sum((data - mean)**2)
    for t in range(0, tmax):
        corr = 0.0
        for i in range(startvalue, len(data) - t):
            corr += (data[i] - mean) * (data[i + t] - mean)
        corr /= norm
        if t > 0:
            int_corr += corr
            sq_sum += corr**2
        corr_err = np.sqrt((1 + 2 * sq_sum) / len(data))

        t_arr.append(t)
        corr_arr.append(corr)
        corr_err_arr.append(corr_err)
        int_corr_err = np.sqrt(sum([i**2 for i in corr_err_arr]))

        int_corr_arr.append(int_corr)
        int_corr_err_arr.append(int_corr_err)

    return t_arr, corr_arr, corr_err_arr, int_corr_arr, int_corr_err_arr


def compute_autocorr_impr_err(data, tmax):
    t_arr, corr_arr, corr_err_arr, int_corr_arr, int_corr_err_arr = corr(data, tmax=tmax, startvalue=0)

    N = len(data)
    err_sqrt = (np.var(data)/N) * (2*int_corr_arr[tmax-1])
    return np.sqrt(err_sqrt)


def compute_autocorr_impr_err_arr(data, tmax):
    I=len(data[0])
    K=len(data[0][0])
    err=np.zeros((I,K))
    mean = np.mean(data, axis=0)
    for i in range(len(data[0])):
        print("Column: ",i)
        for k in range(len(data[0][0])):
            err[i,k]=compute_autocorr_impr_err(data[:,i,k],tmax)
    return np.transpose(mean),np.transpose(err)