# 
# autocorrelation.py                                                               
#
# D. Clarke
# 
# Calculation of integrated autocorrelation time. The main methods have been adapted from software of
# Bernd Berg, Markov Chain Monte Carlo and their Statistical Analysis, World Scientific, 2004, ISBN=978-981-3106-37-6
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import std_mean, std_err, checkTS, binSeries
from latqcdtools.base.plotting import plt, clearPlot, plot_dots, set_params
from latqcdtools.statistics.jackknife import jackknife 
from latqcdtools.base.readWrite import writeTable
from latqcdtools.base.check import checkType


def remove1Jackknife(ts) -> np.ndarray:
    """ 
    Create remove-1 jackknife list from 1-d series.

    Args:
        ts (array-like): time series 

    Returns:
        np.array: 1-d array of jackknife means 
    """
    checkTS(ts)
    return jackknife(np.mean,ts,numb_blocks=len(ts),return_sample=True,conf_axis=0)[0]


def tauint(nt,ts,xhat = None) -> np.ndarray:
    """ 
    Given a time series, calculate estimators for its integrated autocorrelation time  at each Markov time separation.

    INPUT:
         nt--The largest you think tau_int could be.
         ts--Time series array of measurments. Must be taken from equilibrium ensemble so that
             time translation invariance holds. List must be in order of Markov chain generation.
       xhat--True mean of time series (if you know it).

    OUTPUT:
      acint--List of integrated autocorrelation times. 
      """
    checkTS(ts)
    ndat=len(ts)
    if nt>=ndat:
        logger.TBRaise("Need nt<ndat.")
    if xhat is not None:
        x=xhat
    else:
        x=std_mean(ts)
    # Create array of autocovariance
    acov=[]
    for it in range(nt+1):
        numt=ndat-it # number of pairs of time series elements separated by computer time it
        c_it=0.
        for i in range(numt):
            c_it=c_it+(ts[i]-x)*(ts[i+it]-x)
        c_it=c_it/numt
        # Bias correction. This is needed for c(t) to be correct in the limit of uncorrelated data. In principle
        # you don't know ahead of time whether your raw data are effectively correlated. The factor is ndat rather than
        # ndat-it because ndat data points were used to calculate x, which is what matters for the bias.
        if xhat is None:
            c_it=c_it*ndat/(ndat-1.)
        acov.append(c_it)
    acov=np.array(acov)
    # Calculate integrated autocorrelation time. This is equivalent to the following code:
    #    acint=[1.]
    #    for it in range(1,nt+1):
    #        acint.append( acint[it-1] + 2.*acov[it]/acov[0] )
    acint     = np.zeros(nt + 1)
    acint[0]  = 1
    acint[1:] = 1 + 2*np.cumsum(acov[1:]/acov[0])
    return acint


def tauintj(nt,nbins,ts,xhat = None) -> np.ndarray:
    """ 
    Given a time series, calculate jackknife bins of integrated autocorrelation time for each Markov time separation.

    INPUT:
          nt--The largest nt at which you think your estimate for tau_int could lie.
       nbins--The number of jackknife bins.
          ts--Time series array of measurements. Must be taken from equilibrium ensemble so that
              time translation invariance holds. List must be in order of markov chain generation
        xhat--True mean of time series (if you know it).

    OUTPUT:
      acintj--2D list indexed by time, then bin number acintj[it][ibin] 
      """
    checkType('int',nt=nt)
    checkType('int',nbins=nbins)
    checkTS(ts)
    ndat=len(ts)
    if not 1<nbins<ndat:
        logger.TBRaise("Need 1<nbins<ndat. nbins, ndat =",nbins,ndat)
    if xhat is not None:
        x=xhat
    else:
        x=std_mean(ts)
    # For each it, create a list of nbins jackknife measurements of the autocovariance. The measurements are an average
    # over binsize data with separation it. Everything is stored in a 2D array acorj[it][ibin].
    acorj=[]
    for it in range(nt+1):
        numt=ndat-it
        binsize=int(numt/nbins)
        acor=np.zeros(nbins)
        # This is a check that we can't spill into the last bin
        itmax=nbins*binsize-binsize
        if it>=itmax:
            logger.TBRaise("it>=itmax.")
        for ibin in range(nbins):
            i1=ibin*binsize
            i2=(ibin+1)*binsize-1
            acor[ibin] = np.sum( (ts[i1:i2+1] - x)*(ts[i1+it:i2+1+it] - x) )/binsize
        # Bias correction
        if xhat is None:
            acor=acor*ndat/(ndat-1.)
        acorj.append(remove1Jackknife(acor))
    acorj = np.array(acorj)
    # The following code is equivlent to this (but faster since it uses numpy)
    #    acintj=[[1.]*nbins]
    #    for it in range(1,nt+1):
    #        tauintbins=[]
    #        for ibin in range(nbins):
    #            tauintbins.append(acintj[it-1][ibin]+2.*acorj[it][ibin]/acorj[0][ibin])
    #        acintj.append(tauintbins)
    acintj = np.ones((nt+1,nbins))
    increments = 2.*acorj[1:]/acorj[0]
    acintj[1:] = np.cumsum(increments, axis=0) + 1
    return acintj


def getTauInt(ts, nbins, tpickMax, acoutfileName = None, showPlot = False):
    """ 
    Given a time series, return estimates for the integrated autocorrelation time and its error.

    INPUT:
         tpickMax--The largest nt where you think your estimate might become unreliable.
            nbins--The number of jackknife bins (for estimating the error in tau_int)
               ts--Time series array of measurements. Must be taken from equilibrium ensemble so that
                   time translation invariance holds. List must be in order of markov chain generation

    OUTPUT:
          tau_int--Estimate for integrated autocorrelation time.
         tau_inte--Its (jackknife) error bar.
           itpick--The Monte Carlo separation at which this method found its estimate for tau_int. 
    """
    checkTS(ts)

    # Get integrated autocorrelation time list and corresponding jackknife list.
    acint  = tauint (tpickMax,ts,)
    acintj = tauintj(tpickMax,nbins,ts)

    # This block outputs time, tau_int and tau_int error, then gives an estimate of tau_int. When tau_int
    # decreases for the first time, we call that our estimate. This is because we know tau_int should be a 
    # monotonically increasing function of t.
    lmonoton=True
    tau_int=0.
    tau_inte=-1
    itpick=-1
    its = np.arange(tpickMax+1)
    ace = std_err(acintj,axis=1)
    for it in its:
        if lmonoton:
            if not acint[it]<tau_int:
                tau_int=acint[it]
                tau_inte=ace[it]
                itpick=it
            else:  # acint[it] < tau_int ==> tau_int decreased
                lmonoton=False

    if acoutfileName is not None:
        writeTable(acoutfileName,its,acint,ace,header=['t','tau(t)','tau_err(t)'])

    if showPlot:
        clearPlot()
        plot_dots(its,acint,ace,color='black',marker=None)
        set_params(xlabel='conf', ylabel='$\\tau_{\\rm int}$')
        plt.show()
        clearPlot()

    return tau_int, tau_inte, itpick


def gewekeZ(ts, first=0.1, last=0.5, nbins=5) -> float:
    """ 
    Geweke's test of whether a Markov chain has thermalized [1]. If it has, the beginning and the end of the
    time series sample the same distribution, so their means agree. Compare the mean A of the first fraction
    first of ts with the mean B of the last fraction last,
        z = (A - B)/sqrt(err_A^2 + err_B^2).
    For a thermalized chain, z is approximately a standard normal number, so |z| of order 1 is expected, while a
    large |z| means the chain was still drifting at its beginning. Successive measurements are autocorrelated,
    so each error is the error of the mean of nbins bins (batch means). This is valid when a bin is longer than
    the integrated autocorrelation time. Geweke estimated these errors from the spectral density at zero
    frequency instead, which is equivalent. The test checks only the observable in ts, and it cannot notice a
    chain that is stuck in one region the whole time.

    [1] J. Geweke, Evaluating the accuracy of sampling-based approaches to the calculation of posterior moments, 
        in Bayesian Statistics 4, edited by J. M. Bernardo, J. O. Berger, A. P. Dawid, and A. F. M. Smith 
        (Oxford University Press, 1992), pp. 169-193.

    Args:
        ts (array-like): time series of measurements, in the order of the Markov chain
        first (float, optional): fraction of ts at its beginning. Defaults to 0.1.
        last (float, optional): fraction of ts at its end. Defaults to 0.5.
        nbins (int, optional): number of bins in each fraction. Defaults to 5.

    Returns:
        float: z 
    """
    checkTS(ts)
    checkType('int',nbins=nbins)
    if (first <= 0) or (last <= 0) or (first + last > 1):
        logger.TBRaise('Need 0 < first, 0 < last, and first + last <= 1. Got first, last =',first,last)
    n  = len(ts)
    nA = int(first*n)
    nB = int(last*n)
    if (nA < nbins) or (nB < nbins):
        logger.TBRaise('Need at least nbins =',nbins,'measurements in each fraction. Got',nA,'and',nB)
    A    = ts[:nA]
    B    = ts[n-nB:]
    diff = std_mean(A) - std_mean(B)
    err2 = std_err(binSeries(A,nbins))**2 + std_err(binSeries(B,nbins))**2
    if err2 == 0:
        if diff == 0:
            return 0.
        return np.copysign(np.inf, diff)
    return diff/np.sqrt(err2)
