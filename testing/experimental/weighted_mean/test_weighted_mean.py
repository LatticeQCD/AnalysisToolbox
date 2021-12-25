#!/usr/bin/env python3
import numpy as np
from latqcdtools.weighted_average import *
from latqcdtools.statistics import *
from numpy.random import randint, normal


def func(data, weights):
    res = weighted_mean(data[:,0], weights)
    return res


def ref_bootstr(data, err, weights):
    means = []
    for i in range(1000):
        samples = []
        for k in range(len(data)):
            samples.append(normal(data[k], err[k], 1)[0])
        means.append(weighted_mean(samples, weights))
    dev = std_dev(means)
    mean = std_mean(means)
    print("REF", mean, dev)


    
#We can only test by comparing to a non-weighted mean.
data = normal(10, 10, 100)
errors = np.full_like(data, 10)

weights = 1/errors**2

print("weighted_mean, error prop.", weighted_mean(data, weights), np.sqrt(weighted_mean_variance(errors)))
print("weighted mean, unbiased mean error", weighted_mean(data, weights), np.sqrt(unbiased_mean_variance(data, weights)), np.sqrt(unbiased_sample_variance(data, weights)))
print("Standard mean, error", std_mean(data), std_err(data), std_dev(data))
