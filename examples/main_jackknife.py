# 
# main_jackknife.py                                                               
# 
# D. Clarke
# 
# Here we will work out some examples how to use the jackknife method. We will
# generate some fake data and compare different approaches to computing a
# mean and error bar. 
# 

import numpy as np
from latqcdtools.statistics.statistics import error_prop, std_err
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.initialize import TBRNG

# To make random fake data, we need a generator. Here we use what is defined
# as the Toolbox standard RNG, which just renames whatever RNG is the most
# recent recommended one by numpy.
rng = TBRNG()

# Given a time series X and a function, we are interested in unbiased
# estimators of <func(X)> and func(<X>). Feel free to play around with
# different numpy functions
function = np.abs

# IMPORTANT: these wrappers f and g are what gets passed to jackknife.
# The jackknife routine expects a wrapper that computes somewhere a mean.
# Where you put the np.mean depends on what you want to compute

# You should be able to take guidance from what you would write down formally
# with pen and paper. If I want an unbiased estimate of <|x|>, I replace <>
# brackets with np.mean and || with np.abs within the wrapper, leading to
# the following:
def f(dat):
    """ <|x|> """
    return np.mean(function(dat))

# Same game for |<x>|
def g(dat):
    """ |<x>| """
    return function(np.mean(dat))

# We won't look at susceptibilities here, but I include an example susceptibility
# so you can see the versatility of the jackknife
def susc(dat):
    """ <x^2> - <x>^2 """
    return np.mean(dat**2) - np.mean(dat)**2

# Draw from normal distribution with true mean Xhat and sigma Shat
Xhat = 0
Shat = 1
X    = rng.normal(Xhat,Shat,1000)


#
# Now come the jackknife estimators:
#


# Jackknife estimate for <func(x)>
m, e = jackknife(f,X)
print('<func(X)> jack',get_err_str(m,e))

# This is not only a naive estimator for <func(x)> but also
# a severely biased estimator for func(<X>), i.e. wrong
m = f(X)
e = std_err(function(X))
print('<func(X)> naiv',get_err_str(m,e))

# Jackknife estimate for |func(X)|
m, e = jackknife(g,X)
print('func(<X>) jack',get_err_str(m,e))

# Naive estimate for func(<X>) using error propagation. 
m, e = error_prop(function,np.array([np.mean(X)]),np.array([std_err(X)]))
print('func(<X>) naiv',get_err_str(m,e))
print('func(<X>) true',function(Xhat))


# This is the most basic usage of jackknife, but one can control the number
# of jackknife bins (numb_blocks). Also the wrappers f and g can be arbitrarily
# complicated: They don't have to just wrap simple numpy operations. If the
# operation is sufficiently time-intensive, you can parallelize the jackknife
# by adjusting the number of processors (nproc). I recommend using nproc=1 almost
# always: the compute-time demand of the wrapper has to be large enough to offset
# the overhead from using multiple processors to justify nproc>1.

