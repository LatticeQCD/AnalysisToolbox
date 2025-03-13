# 
# main_isingModel.py 
# 
# D. Clarke 
# 
# In this example, we show how one can use the AnalysisToolbox to carry out a general statistical physics
# computation. We take as our exemplar the Ising model in N dimensions, which is more or less the
# most simple and renowned statistical mechanical model. We simulate it using Markov chain Monte Carlo.
#

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


initialize()
latexify()


#
# Simulation parameters. 
#
Nd    = 2        # number of dimensions
Tlow  = 2.25     # lowest temperature to sample (kB=1)
Thi   = 2.4      # highest temperature to sample
h     = 0.       # external magnetic field
L     = 8        # spatial extension
Nequi = 2000     # equilibrate with this many MCMC sweeps
Nmeas = 600      # measure this many MCMC sweeps
Nskip = 5        # separate measurements by this many MCMC sweeps
start = 'hot'    # start with all spins up (up), down (down), or random (hot)


Tlist = np.linspace(Tlow,Thi,DEFAULTTHREADS)
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
    logger.info('T =',T)

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


    def actionLocal(coord):
        """ Local contribution to action.

        Args:
            coord (array-like): local coordinate 
        """
        s0 = lat.getElement(coord)
        NNsum = 0.
        for mu in range(Nd):
            NNsum += lat.getElement(lat.march(coord,mu,1)) + lat.getElement(lat.march(coord,mu,-1))
        return s0*NNsum


    def getAction():
        return lat.bulkReduce(actionLocal)
        

    def mcLocal(coord):
        """ The Markov step. 

        Args:
            coord (array-like): local coordinate 
        """
        s0 = lat.getElement(coord)
        dH = 2*beta*actionLocal(coord)
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
    actions        = []
    for imeas in range(Nmeas):
        for iskip in range(Nskip):
            lat.iterateOverRandom(mcLocal)
        magnetizations.append(getMagnetization())
        actions.append(getAction())
    magnetizations = np.array(magnetizations) 
    actions        = np.array(actions) 


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

    writeTable(f'ising_L{L}_d{Nd}/T{T}.d',magnetizations,actions)

    return Mm, Me, chim, chie, Bm, Be


# Parallelize over the temperatures. This parallelization strategy is easiest since each temperature
# amounts to an independent run.
data = parallel_function_eval(runIsingModel,Tlist)


res_M, res_E, res_chi, res_chie, res_B, res_Be = [], [], [], [], [], []
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
           header=['T','|M|','|M|_err','chi','chi_err','B','B_err'])
finalize()


