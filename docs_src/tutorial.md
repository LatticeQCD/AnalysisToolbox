# Tutorial

Here we walk through some examples found in the `latqcdtools/examples`
directory. We try to showcase the flexibility of things like the
[bootstrap](dataAnalysis/bootstrap.md) routine, some convenience wrappers
for plotting, and some pre-packaged physics analysis code.

## Hadron resonance gas calculation

The hadron resonance gas (HRG) model, and its implementation in the AnalysisToolbox,
is described in some detail [here](physicsAnalysis/HRG.md). Below is a small code
that highlights some of the features of the AnalysisToolbox. For example
`initialize` simulataneously saves all screen output to the log file `HRG.log`
along with the git commit hash, so you can track down which version of the
AnalysisToolbox you used.

Here is `latqcdtools/examples/main_HRG_simple.py`

```Python
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.HRG import HRG
from latqcdtools.base.readWrite import readTable, writeTable
from latqcdtools.base.initialize import initialize, finalize

# Write terminal output to log file. Includes git commit hash.
initialize('HRG.log')

# Pick a temperature range in MeV
T = np.arange(100, 166, 1)

# Read in hadron names, masses, charges, baryon number, strangeness,
# charm, and degeneracy factor. This table is provided with AnalysisToolbox.
QMHRG_table = '../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt'
hadrons, M, Q, B, S, C, g = readTable(QMHRG_table, usecols=(0,1,2,3,4,5,6),
                                      dtype="U11,f8,i8,i8,i8,i8,i8")
w = np.array([1 if ba==0 else -1 for ba in B])

# Instantiate HRG object.
QMhrg = HRG(M,g,w,B,S,Q,C)

# This computation is vectorized since T is a numpy array.
logger.info('Computing chi2B.')
chi = QMhrg.gen_chi(T, B_order=2, Q_order=0, S_order=0, C_order=0,
                    muB_div_T=0.3, muQ_div_T=0, muS_div_T=0, muC_div_T=0)

# Output T and chi2B as columns in this table.
writeTable("chi2B.txt", T, chi, header=['T [MeV]','chi2B (QMHRG)'])

finalize()
```

In this case we work at fixed $\mu_B/T=0.3$
and compute $\chi_2^B$ as function of temperature. What is needed is as much
relevant input knowledge about hadron bound states as is known; this is
collected by `readTable`. This is a wrapper for `numpy.loadtxt()`, hence you
see you can pass it many of the same keyword arguments. The next line
for each species whether the gas is bosonic or fermionic.

Finally the `HRG` class is instantiated as `QMhrg`, which besides
generic conserved charge cumlants like the one computed with `gen_chi`,
contains many methods for various thermodynamic observables such as the
pressure and entropy.
The results are saved in a table `chi2B.txt`.

## Ising model

The AnalysisToolbox is equipped with a `Lattice` class. This class takes care of
indexing, assumes you have periodic boundary conditions, and has iterators, i.e.
functions that perform some operation on every site of a subset of the lattice.
This allows you to write your own statistical mechanical simulations. 

Here is `latqcdtools/examples/main_isingModel.py`

```Python

import numpy as np
from latqcdtools.physics.lattice import Lattice
from latqcdtools.base.initialize import initialize, finalize
from latqcdtools.base.plotting import latexify, plt, plot_dots, set_params, clearPlot
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.statistics.statistics import std_mean
from latqcdtools.base.speedify import parallel_function_eval, DEFAULTTHREADS
from latqcdtools.base.readWrite import writeTable
from latqcdtools.statistics.jackknife import jackknife
import latqcdtools.base.logger as logger


initialize('example_ising.log')
latexify()


#
# Simulation parameters. 
#
Nd    = 2        # number of dimensions
Tlow  = 1.0      # lowest temperature to sample (kB=1)
Thi   = 3.0      # highest temperature to sample
h     = 0.       # external magnetic field
L     = 8        # spatial extension
Nequi = 300      # equilibrate with this many MCMC sweeps
Nmeas = 100      # measure this many MCMC sweeps
Nskip = 5        # separate measurements by this many MCMC sweeps
start = 'hot'    # start with all spins up (up), down (down), or random (hot)


Tlist = np.linspace(Tlow,Thi,2*DEFAULTTHREADS)
logger.info()
logger.info('Starting parameters:')
logger.info('  Nequi =',Nequi)
logger.info('  Nmeas =',Nmeas)
logger.info('  Nskip =',Nskip)
logger.info('   Latt =',str(L)+'^'+str(Nd))
logger.info('      h =',h)
logger.info()


#
# We wrap our simulation inside a function. This allows us to trivially parallelize our run
# over the temperatures using parallel_function_eval.
#
def runIsingModel(T):

    beta  = 1/T
 
    # Initialize the random number generator. When carrying out a statistical physic MCMC, it's
    # crucially important that you pick a good one. The default_rng() constructor is what
    # numpy recommends, which at the time of writing utilizes O'Neill's PCG algorithm. 
    rng = np.random.default_rng()

    # Initialize the lattice object. The first argument says the lattice geometry will be L**Nd,
    # while the second argument is an example of the kind of object that will be on each site.
    # Here each site has a scalar, so I just put in the number 1.
    lat = Lattice( (L,)*Nd, 1 )

    def initialize_site(coord):
        if start == 'up':
            lat.setElement(coord,1)
        elif start == 'down':
            lat.setElement(coord,-1)
        elif start == 'hot':
            lat.setElement(coord,rng.choice([1,-1]))
        else:
            logger.TBError('Unknown start option',start)

    # We then initialize all sites. Here interateOverRandom, whose naming convention follows that
    # of SIMULATeQCD, is an iterator that applies the function initialize_site to each site.
    # The function should take an array-like coord as argument. Look into the Lattice class to
    # see what other iterators are possible. 
    lat.iterateOverRandom(initialize_site)
    

    def getMagnetization():
        """ Order parameter is |M|. """
        return np.abs(np.sum(lat.grid)/lat.vol)
    
    
    def mcLocal(coord):
        """ The Markov step. 

        Args:
            coord (array-like): local coordinate 
        """
        s0 = lat.getElement(coord)
        NNsum = 0.
        for mu in range(Nd):
            NNsum += lat.getElement(lat.march(coord,mu,1)) + lat.getElement(lat.march(coord,mu,-1))
        dH = 2*beta*s0*NNsum
        if dH < 0:
            lat.setElement(coord,-s0)
        elif rng.random() < np.exp(-dH):
            lat.setElement(coord,-s0)
    

    # Equilibrate the lattice by carrying out Nequi equilibration steps.   
    for iequi in range(Nequi):
        lat.iterateOverRandom(mcLocal)
    

    # With the remaining MCMC steps, we measure the magnetization. We skip every Nskip MCMC 
    # steps to reduce autocorrelation.    
    magnetizations = []
    for imeas in range(Nmeas):
        for iskip in range(Nskip):
            lat.iterateOverRandom(mcLocal)
        magnetizations.append(getMagnetization())
    magnetizations = np.array(magnetizations) 


    def chi_M(M):
        """ Magnetic susceptibility. """
        return lat.vol*( std_mean(M**2)-std_mean(M)**2 )


    def binder(M):
        """ Binder cumulant. For the ising model in 2d, this should tend to 2/3 below
        Tc, while it should tend to 0 above Tc, in the thermodynamic limit. """
        return 1-np.mean(M**4)/(3*np.mean(M**2)**2)


    # We have parallelized over the temperatures, so we're not allowed to parallelize over the jackknife:
    # that would be nested parallelization!
    Mm  , Me   = jackknife( std_mean, magnetizations, nproc=1 )
    chim, chie = jackknife( chi_M   , magnetizations, nproc=1 )
    Bm  , Be   = jackknife( binder  , magnetizations, nproc=1 )
    logger.info('T, <|M|>, chi, B =',round(T,2),get_err_str(Mm,Me),get_err_str(chim,chie),get_err_str(Bm,Be))
    return Mm, Me, chim, chie, Bm, Be


# Parallelize over the temperatures. This parallelization strategy is easiest since each temperature
# amounts to an independent run.
data = parallel_function_eval(runIsingModel,Tlist)


res_M    = []
res_E    = []
res_chi  = []
res_chie = []
res_B    = []
res_Be   = []
for i in range(len(Tlist)):
    res_M.append(   data[i][0])
    res_E.append(   data[i][1])
    res_chi.append( data[i][2])
    res_chie.append(data[i][3])
    res_B.append(   data[i][4])
    res_Be.append(  data[i][5])


# Plot the magnetization.
plot_dots(Tlist,res_M,res_E)
set_params(xlabel='$T$',ylabel='$\\ev{|M|}$',title='$V='+str(L)+'^'+str(Nd)+'$ Ising model',alpha_xlabel=0)
plt.show()
clearPlot()


# Plot the magnetic susceptibility.
plot_dots(Tlist,res_chi,res_chie)
set_params(xlabel='$T$',ylabel='$\\ev{\\chi}$',title='$V='+str(L)+'^'+str(Nd)+'$ Ising model',alpha_xlabel=0)
plt.show()


# Record the results in a nicely formatted table.
writeTable('ising_'+str(L)+'_'+str(Nd)+'.d',Tlist,res_M,res_E,res_chi,res_chie,res_B,res_Be,
           header=['T','|M|','|M|_err','chi','chi_err','B','Be'])
finalize()
```

Again, the most useful feature of this example is the `Lattice` class. This wraps a
lattice that saves an object of arbitrary type at each site, here scalars. You can move
one site in a direction using `march`. In the Metropolis step, we feature the lattice
accessors `getElement` and `setElement`. To deal with autocorrelation and propagate
error, we use the [jackknife](dataAnalysis/jackknife.md) class.

## Continuum-limit extrapolation

Suppose we want to perform a continuum-limit extrapolation to
determine the deconfinement transition temperature $T_d$ in pure SU(3). The order
parameter for this phase transition is given by the Polyakov loop, $P$. 
The transition is first-order in the
thermodynamic limit, where $\langle |P|\rangle$ as function of temperature
would jump discontinuously at $T_d$.
At finite volume, this abrupt jump becomes smooth, and $T_d$ is estimated by the
inflection point of the curve.

Here is `latqcdtools/examples/main_continuumExtrapolate.py`

```Python
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.initialize import initialize, finalize
from latqcdtools.math.num_deriv import diff_deriv
from latqcdtools.math.spline import getSpline
from latqcdtools.statistics.statistics import gaudif
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.physics.continuumExtrap import continuumExtrapolate
from latqcdtools.physics.constants import r0_phys
from latqcdtools.physics.lattice_params import latticeParams

initialize('cont.log')

Nts         = [6,8,10,12,14,16,18,20]
Tds, Tderrs = [], []

for Nt in Nts:

    T  = []
    Ns = Nt*3

    # Read in Polyakov loop measurements, 
    beta, PM, PE = readTable('ploop/Nt'+str(Nt)+'.txt',usecols=(0,1,2))

    # Create array of temperatures in physical units
    for b in beta:
        lp = latticeParams(Ns, Nt, b, scaleType='r0') 
        T.append( lp.getT() )
    t = np.linspace(T[0],T[-1],1001)

    # Extract Td from inflection point of <|P|> vs T using natural spline
    def getTd(pm):
        spl  = getSpline(T, pm, natural=True)
        dPdT = diff_deriv(t, spl)
        maxIndex = np.argmax(dPdT)
        return t[maxIndex]

    # Error in Td estimate comes from 1000 Gaussian bootstrap samples
    Td, Tde = bootstr_from_gauss(getTd, PM, PE, 1000)
    Tds.append(Td)
    Tderrs.append(Tde)

# Perform O(a^4) continuum-limit extrapolation
result, result_err, chidof = continuumExtrapolate( Nts, Tds, Tderrs, order=2, xtype="Nt",
                                                   show_results=True, plot_results=True )

# Do a Z-test against literature result,  
r0 = r0_phys(year=2014, units="MeVinv")
Tdr0, Tdr0e = r0 * result[0], r0 * result_err[0]
Tdr0_lit, Tdr0_lite = 0.7457, 0.0045
logger.info('q(ours vs. lit) =',gaudif(Tdr0,Tdr0e,Tdr0_lit,Tdr0_lite))

finalize()
```

Above we show how such an extrapolation is achieved with
the AnalysisToolbox, along with error estimation, plotting the results, and
carrying out a statistical comparison with the known literature value.
We assume you already have results for $\langle |P|\rangle$ at various $N_\tau$, which we
read in from tables of the form `Nt6.txt`. 
For each $N_\tau$, this code estimates the inflection point of $\langle |P|\rangle$ as a function 
of $T$ by taking a [derivative](math/calculus.md) using a 
[spline](dataAnalysis/curveFitting.md) to get $T_d(N_\tau)$. 
Temperatures are calculated in MeV with the help of our
[class](physicsAnalysis/latticeParameters.md) that collects ensemble parameters. This 
procedure is wrapped in a user-defined function `getTc`, 
so that errors in the $\langle |P| \rangle$ data can be 
conveniently propagated into the error in $T_d(N_\tau)$
using a Gaussian [bootstrap](dataAnalysis/bootstrap.md). 

Having the `Nts`, `Tds`, and `Tderrs`, we are ready to perform a
continuum-limit extrapolation. This will perform an
extrapolation to second order in $a^2$, i.e. $\mathcal{O}(a^4)$, print the fit
results to screen, and create a plot of the extrapolation for you.
The arrays `result` and `result_err` contain the best fit parameters
along with their errors, with `result[0]` being the continuum value $T_d$.
In the last block we compare our result with $T_d r_0$ from the
[literature](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.91.096002). 
The temperatures calculated in this code implicitly had units of
MeV, hence we need $r_0$ in [physical units](physicsAnalysis/referenceScales.md). 
Finally we call `gaudif` to carry out a Gaussian difference test or
Z-test, which is implemented in our [statistics](dataAnalysis/statistics.md) module.

## Statistical bootstrap

Here is `latqcdtools/examples/main_bootstrap.py`

```Python
import numpy as np
from latqcdtools.base.readWrite import readTable
from latqcdtools.math.spline import getSpline
from latqcdtools.base.plotting import plt, plot_dots, plot_band, latexify
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.base.initialize import initialize, finalize

initialize('example_bootstrap.log')

latexify()

xdata, ydata, edata = readTable("../testing/statistics/wurf.dat", usecols=(0,2,3))

plot_dots(xdata,ydata,edata)

# We will interpolate to these x-values
xspline = np.linspace(np.min(xdata),np.max(xdata),101)

def splineError(data):
    """ We assume no errors in the xdata. For the ydata, we will pass the bootstrap
    routine ydata along with the errors. The input to this function, data, will then
    be generated in each bootstrap bin by drawing normally from ydata with a spread
    of edata. In this example we simply get a spline, but you wrap anything you want
    inside your bootstrap procedure.
    """
    ys = getSpline(xdata,data,3)
    return ys(xspline)

ybs, ybserr = bootstr_from_gauss(splineError,data=ydata,data_std_dev=edata,numb_samples=100)

plot_band(xspline,ybs-ybserr,ybs+ybserr,label='bootstrapped interpolation')

plt.show()

finalize()
```

In the above example we created an interpolated band using the data; there are 100 bootstrap
samples, where each data point is drawn from `normal(ydata,edata)`. The error is taken by
default to be the 68-percentile bounds, and the central value is given as the median. Like
the `jackknife` from the Ising model example, the bootstrap routine is parallelized by default.
