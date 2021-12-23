#!/usr/bin/env python3
import numpy as np
from numpy.random import normal
from latqcdtools.tools import *
from latqcdtools.statistics import *
import latqcdtools.bootstr as bs
import latqcdtools.jackknife as jk


data = normal(2, 1, (1000, 2))

data[:,1] = normal(4, 2, 1000)

def ref_std_dev(data, axis):
    mean = np.mean(data, axis = axis)
    error = np.zeros_like(mean)
    data = np.rollaxis(data, 0, axis + 1)

    for i in range(0, len(data)):
        error += (data[i] - mean) * (data[i] - mean)
    error = np.sqrt(error / (len(data) - 1))
    return error


def ref_std_err(data, axis):
    data = np.rollaxis(data, 0, axis + 1)
    return std_dev(data, 0) / np.sqrt(len(data))



print_results(std_dev(data, 0), ref_std_dev(data, 0), text = "Standard deviation, axis = 0")

data = normal(2, 1, (2, 1000))
data[1:] = normal(4, 2, 1000)

print_results(std_dev(data, 1), ref_std_dev(data, 1), text = "Standard deviation, axis = 1")

print_results(std_err(data, 1), ref_std_err(data, 1), text = "Standard error, axis = 1")



print("Testing bootstrap and jackknife. Please compare the following results")

print("Reference from std_mean and std_err")
print(mean_and_err(data, axis = 1))

print("\nBootstrap with tuples")
mean, err = bs.bootstr(lambda x: (np.mean(x[0]), np.mean(x[1])), data, 1000)
print(mean, err)

print("\nBootstrap with tuples, error by distribution")
mean, err = bs.bootstr(lambda x: (np.mean(x[0]), np.mean(x[1])), data, 1000, 
        err_by_dist = True)
print(mean, err)

def tupleTester(x):
    return (np.mean(x[0]),np.mean(x[1]))

print("\nJackknife with tuples")
mean, err = jk.jackknife(tupleTester, data)
print(mean, err)

print("\nBootstrap with array")
mean, err = bs.bootstr(np.mean, data, 1000, args = {'axis' : 1})
print(mean, err)

print("\nJackkife with array")
mean, err = jk.jackknife(np.mean, data, args = {'axis' : 1})
print(mean, err)