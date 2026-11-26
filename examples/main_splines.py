# 
# main_splines.py                                                               
# 
# D. Clarke
# 
# Let's go ahead and see how splines can be used in the toolbox. We're going to look at
# two data sets. The first plots the magnetic susceptibility against the temperature in
# an O(3) model. The second plots the QCD trace anomaly against the temperature in
# Nf=1+1+1 QCD. 
# 

import numpy as np

from latqcdtools.base.plotting import plot_dots, set_params, plt, clearPlot, plot_lines, \
    plot_fill, plot_vspan
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.speedify import DEFAULTTHREADS
import latqcdtools.base.logger as logger

from latqcdtools.math.spline import getSpline, Interpolating, Regression, KnotAverage, \
    Smoothing, bootSpline, posteriorSpline
from latqcdtools.math.num_deriv import diff_deriv

#
# First let's look at the O(3) model data. Notice there is sharp curvature in the middle,
# which flags the pseudocritical temperature Tc(L=48). More info about how to use the
# plotting class is in main_plotting.py.
#
T, _, _, chi, chie = readTable('../datasets/example_files/heisenberg3d_L48_chi.txt')
plot_dots(T,chi,chie,marker=None)
set_params(xlabel='$T$',ylabel='$\\chi_{|m|}$')
plt.show()

#
# There are several different kinds of splines available in the Toolbox:
#
#   Interpolating--Natural spline with interior knots fixed at each datum. This
#                  spline passes through each datum.
#      Regression--You specify the knots ahead of time. These can be evenly spaced in
#                  the data index, randomly chosen, specified by the user, or a combination
#                  of user-specified and evenly-spaced/random. In this latter case, the
#                  automatic choices are made BEFORE looking at the user-specified knots.
#                  Then, the regression spline minimizes chi^2.
#       Smoothing--Interior knots at each datum, and a smoothing penalty is included. In
#                  particular, splines with high global curvature are punished.
#     KnotAverage--Here, a Markov chain of possible numbers of knots and knot locations
#                  is created. The accept/reject is determined through the AIC, i.e.
#                  AIC/2 plays the role of an action. The final result is an average
#                  over the visited splines.
#
# Under the hood, every spline in the Toolbox is eventually represented in the B-spline
# basis as a TBSpline class. This allows a unified way of handling splines, no matter
# how the spline was arrived at. Also under the hood, every spline is hard-coded to be
# cubic, except for regression splines. Regression splines are cubic by default, but the
# user has the option of changing the degree.
#

#
# ----------------------------------------------------------------------- BASIC SPLINES
# 

#
# With this discussion out of the way, let's try some naive spline fits to the data. If
# we aren't interested in uncertainties, getSpline() is the most convenient method. It
# returns the spline representation of the data.
# 
plot_dots(T,chi,chie,marker=None)
set_params(xlabel='$T$',ylabel='$\\chi_{|m|}$')

spl = getSpline(xdata=T, ydata=chi, kind=Interpolating())

#
# A spline representation spl from getSpline can be called like an ordinary python
# function, i.e. spl(x) evaluates the spline at x. This makes it convenient to plug
# into other Python functions, like plotting functions 
#
Tgrid = np.linspace(np.min(T), np.max(T), 1001)
plot_lines(Tgrid,spl(Tgrid),marker=None,label='interp.')
plt.show()

clearPlot()

#
# If you want, you can also look at its nth derivative, spl.n_deriv(Tgrid,n). 
# 
# Now let's compare some different kinds of regression splines.
#
plot_dots(T,chi,chie,marker=None)
set_params(xlabel='$T$',ylabel='$\\chi_{|m|}$')

# Cubic, knots evenly spaced in data index
spl = getSpline(xdata=T, ydata=chi, edata=chie, kind=Regression(nknots=8))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='reg. even-spaced')

# Cubic, random knots 
spl = getSpline(xdata=T, ydata=chi, edata=chie, kind=Regression(nknots=8,rand=True,seed=1))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='reg. random')

# Cubic, I choose two knots, and the other 6 are evenly spaced in data
spl = getSpline(xdata=T, ydata=chi, edata=chie, kind=Regression(nknots=8,knots=[1.45,1.5]))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='reg. my knots')

# Maybe try a quartic one 
spl = getSpline(xdata=T, ydata=chi, edata=chie, kind=Regression(nknots=8,order=4))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='reg. quartic even-spaced')
plt.show()

clearPlot()

#
# As you can see, these all perform rather poorly. They could be made to work
# with a lot of trial-and-error, through some combination of changing the
# interpolation range, number of knots, spline degree, and knot location, but
# that's not really ideal.
#
# We can take out the guesswork with the knots at least by performing a
# Monte Carlo over knot locations. Let's see how that works. We can optionally
# specify:
#       kmax--Try at most 30 knots
#      nequi--Perform 10000 equilibration/thermalization updates
#   nupdates--Average over 4000 updates after equilibration 
#
plot_dots(T,chi,chie,marker=None)
set_params(xlabel='$T$',ylabel='$\\chi_{|m|}$')

splMC = getSpline(T, chi, edata=chie, kind=KnotAverage(kmax=30,nequi=10000,nupdates=4000))
plot_lines(Tgrid,splMC(Tgrid),marker=None,label='MCMC')

# For comparison, we show also the interpolating spline, and create an inset
# that highlights the differences between these methods
splInterp = getSpline(xdata=T, ydata=chi, kind=Interpolating())
plot_lines(Tgrid,splInterp(Tgrid),marker=None,label='interp.')

inset = plt.axes([0.18, 0.42, 0.26, 0.38])
plot_dots(T,chi,chie,marker=None,ax=inset)
plot_lines(Tgrid,splMC(Tgrid),marker=None,ax=inset)
plot_lines(Tgrid,splInterp(Tgrid),marker=None,ax=inset)
set_params(xmin=1.26,xmax=1.385,ymin=2.4,ymax=4.9,ax=inset)
plt.show()

clearPlot()

#
# We see the MCMC spline does quite well, and the inset shows us that the
# MCMC spline is not really overfitting! Finally, let's have a look at
# a smoothing spline. The parameter lambda weights the importance of the
# curvature punishment. We plot for a few different possible lambda.
#
plot_dots(T,chi,chie,marker=None)
set_params(xlabel='$T$',ylabel='$\\chi_{|m|}$')

spl = getSpline(T, chi, edata=chie, kind=Smoothing(lam=0.1))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='$\\lambda=0.1$')

spl = getSpline(T, chi, edata=chie, kind=Smoothing(lam=0.01))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='$\\lambda=0.01$')

spl = getSpline(T, chi, edata=chie, kind=Smoothing(lam=0.001))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='$\\lambda=0.001$')

spl = getSpline(T, chi, edata=chie, kind=Smoothing(lam=0.0001))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='$\\lambda=0.0001$')

spl = getSpline(T, chi, edata=chie, kind=Smoothing(lam=1e-12))
plot_lines(Tgrid,spl(Tgrid),marker=None,label='$\\lambda\\approx0$')
plt.show()

clearPlot()

#
# As lambda-->0, the spline approaches an interpolating spline. Since
# this spline punishes global curvature, it is generally not a good
# choice for data like this susceptibility, which exhibits a large
# peak locally in the data. Again, this could be solved by adjusting
# the interpolation range, but when making splines, we would like
# to avoid as many trial-and-error "model choices" as we can.
#

#
# ---------------------------------------------------------- PROPAGATING UNCERTAINTIES 
#

#
# getSpline() is our basic spline method, but it doesn't provide a way to compute
# uncertainties by itself. The Toolbox has two other wrappers that are useful for
# propagating the uncertainty in the data to uncertainty in the spline. We will
# briefly look at both. 
#
# A common strategy for propagating uncertainty when original measurements aren't
# available is the bootstrap. The bootSpline() method uses a parametric (Gaussian)
# bootstrap to propagate uncertainties from the data. This is justified whenever
# the data were estimated from the standard arithmetic mean due to the central
# limit theorem.
#
# The above examples can therefore be easily extended by swapping out the getSpline()
# wrapper with bootSpline(). This method requires the number of points on the x-grid,
# nsupport. Now, what is returned is a dictionary res that contains
# a lot of useful information:
#     res['splineMean']--Analogue of spl above; does not contain error information
#     res['xspl']      --support points
#     res['yspl']      --bias-corrected spl on the support points xspl
#     res['ysple']     --error in spl evaluated on the support points xspl
#     res['xmax']      --x-value at which y takes its maximum 
#     res['xmaxe']     --uncertainty in xmax
# You can look at the bootSpline() docstring to see other possibilities. Here we
# showcase a bootSpline() of the KnotAverage kind, which will, in each bootstrap
# sample, carry out the spline Monte Carlo.
#
plot_dots(T,chi,chie,marker=None)
set_params(xlabel='$T$',ylabel='$\\chi_{|m|}$')

# Doing MCMC inside a bootstrap loop is generally a bit slow. Here we speed it up
# by parallelizing the bootstrap procedure. DEFAULTTHREADS is an automatic choice
# based on the number of threads available on your computer. Here, we use
# Nboot=200 bootstrap samples. It takes a few minutes on my laptop.
logger.info('Doing bootstrap+MCMC...')
res = bootSpline(T,chi,chie,kind=KnotAverage(kmax=30,nequi=10000,nupdates=4000),
                 numb_samples=200,nsupport=1001,nproc=DEFAULTTHREADS)

# We plot the spline with its uncertainty using the res dictionary
plot_fill(res['xspl'],res['yspl'],res['ysple'])

# Also plot the peak location with its uncertainty. In the context of the 3-d O(3)
# model, this serves as an estimate of the pseudocritical temperature Tc(L=48)
Tc  = res['xmax']
Tce = res['xmaxe']
plot_vspan(Tc-Tce,Tc+Tce,color='red',alpha=0.3)
logger.info(f'Tc(L=48) = {get_err_str(res["xmax"],res["xmaxe"])}')
plt.show()

clearPlot()

#
# Great, that worked well. I invite you to zoom in to different regions of the
# plot to see how this spline performed.
#
# Finally, let's look at a different data set. Here, we want to do a spline fit
# to the trace anomaly in QCD. Such spline fits are needed to compute the QCD
# pressure, which is more or less the integral over the trace anomaly. In this
# context, one needs to know the value of the trace anomaly at some low
# temperature. Here, we take this information from the HRG model. (You can learn
# more about the HRG in main_HRG_simple.py.) 
#
set_params(xlabel='$T$ [MeV]',ylabel='$\\Theta^{\\mu\\mu}/T^4$')

# Set up our HRG model
from latqcdtools.physics.HRG import HRG
_, M, Q, B, S, C, g = readTable('../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt',
                                usecols=(0,1,2,3,4,5,6), dtype="U11,f8,i8,i8,i8,i8,i8")
hrg = HRG(M,g,np.where(B==0,1,-1),B,S,Q,C)

# Our low-temperature anchor point is at 100 MeV. 
T0 = 100

# The trace anomaly arrays will be the measured data with the 100 MeV reference
# value inserted as the first element of each array. The trace anomaly computed
# in the HRG will be assigned an artificial 10% uncertainty.
T, TA, TAe = readTable('../datasets/example_files/traceAnomaly.txt')
T   = np.insert(T  ,0,                     T0 )
TA  = np.insert(TA ,0,    hrg.theta_div_T4(T0))
TAe = np.insert(TAe,0,0.1*hrg.theta_div_T4(T0))
plot_dots(T,TA,TAe,marker=None)

# Besides using this HRG datum, we will also ask the spline to have its derivative
# at 100 MeV agree with the slope of the trace anomaly there. The first derivative
# can be constrained for these splines by passing a deriv argument 
slope  = diff_deriv(T0,hrg.theta_div_T4)

# Here I show two splines: the bootstrap-MCMC from earlier, where nequi and nupdates
# take default values, and a "posterior spline". The posterior spline is a smoothing
# spline that averages over many possible smoothing parameters lambda, i.e. the
# lambda is "integrated out". Means and uncertainties for the posterior spline are 
# computed using Bayesian model averaging, so the uncertainty includes both statistical
# uncertainty along with a systematic part that measures the spread due to different
# choices of lambda.
resMC  = bootSpline(     T, TA, TAe, kind=KnotAverage(),nsupport=101 ,deriv=(T0,slope,0.1*slope),nproc=DEFAULTTHREADS)
resPOS = posteriorSpline(T, TA, TAe, kind=Smoothing() , nsupport=101, deriv=(T0,slope,0.1*slope))
plot_fill(resMC['xspl'] ,resMC['yspl'] ,resMC['ysple'] ,color='green', label='MCMC')
plot_fill(resPOS['xspl'],resPOS['yspl'],resPOS['ysple'],color='red'  , label='posterior')
plt.show()

#
# The two splines more or less overlap, except at the highest two temperatures. There,
# the Monte Carlo tries to pass through the 163 MeV point. This leads to extra curvature
# between the highest two temperatures, and the MCMC tries to put a maximum there, even
# though there is no evidence from the data that there should be such a maximum; indeed
# the trace anomaly at 174 MeV is larger than at 163 MeV. The posterior spline handles
# this situation better due to its curvature penalty, which makes it harder for splines
# to add additional structure between data.
#
# To summarize: when choosing between these two types of splines, I generally recommend 
# the posterior spline, which is faster, avoids extra structure, and is model-independent
# to the extent that it doesn't depend on lambda. The exception is the failure mode we saw
# with smoothing splines on the 3-d O(3) data: it doesn't do well when the data have a large, local curvature,
# like what a susceptibility has near a phase transition.
#
