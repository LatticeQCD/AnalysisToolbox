import numpy as np
from numpy.linalg import inv
from latqcdtools.statistics import std_median, dev_by_dist
from numpy.random import normal
import math
    

def cov_weighted_average(data, cov):
    """Compute the weighted mean if correlations are present."""
    W = np.ones_like(data)
    sig2 = 1 / (W.transpose().dot(inv(cov).dot(W)))
    av = sig2 * (W.transpose().dot(inv(cov).dot(data)))
    return av, np.sqrt(sig2)



def chi_square(data, weights = None):
    if weights is None:
        weights = np.ones_like(data)
    data = np.asarray(data)
    mean = weighted_mean(data, weights)
    return np.sum(weights.dot((data - mean)**2))



def bootstr_add_dist(data, errors, nstat = 1000, plot_hist = False):
    """Given possibly correlated data and corresponding error bars that are assumed to be Gaussian,
    resample as follows: For data point i, draw nstat new resampled measurements from a
    Gaussian distribution with mean data[i] and standard deviation errors[i]. Concatenate
    these ndata * nstat resampled measurements into a new distribution. Return the median of
    this new distribution. The error is taken as the distance between the median and the
    68% percentile. Return that, too."""
    dists = ()
    for i in range(len(data)):
        dist = normal(data[i], errors[i], size = nstat)
        dists += (dist,)
    tot_dist = np.concatenate(dists)
    med, err = std_median(tot_dist), dev_by_dist(tot_dist)
    err_dist = tot_dist[(tot_dist < med + err) & (tot_dist > med - err)]
    if plot_hist:
        hist, bins = np.histogram(tot_dist, bins = math.ceil(nstat / 10))
        plt.hist(tot_dist, bins = bins)
        plt.hist(err_dist, bins = bins)
    return med, err
