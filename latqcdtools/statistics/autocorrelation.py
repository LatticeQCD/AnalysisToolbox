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
from latqcdtools.statistics.statistics import std_mean, std_err, checkTS
from latqcdtools.base.plotting import plt, clearPlot, plot_file
from latqcdtools.statistics.jackknife import jackknife 


def remove1Jackknife(ts):
    """ Create remove-1 jackknife list from 1-d series.

    Args:
        ts (array-like): time series 

    Returns:
        np.array: 1-d array of jackknife means 
    """
    checkTS(ts)
    return jackknife(np.mean,ts,numb_blocks=len(ts),return_sample=True,conf_axis=0)[0]


def tauint(nt,ts,xhat = None):
    """ Given a time series, calculate estimators for its integrated autocorrelation time  at each Markov time separation.

    INPUT:
         nt--The largest you think tau_int could be.
         ts--Time series array of measurments. Must be taken from equilibrium ensemble so that
             time translation invariance holds. List must be in order of Markov chain generation.
       xhat--True mean of time series (if you know it).

    OUTPUT:
      acint--List of integrated autocorrelation times. """
    checkTS(ts)
    ndat=len(ts)
    if nt>=ndat:
        logger.TBError("Need nt<ndat.")
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
    # Calculate integrated autocorrelation time
    acint=[1.]
    for it in range(1,nt+1):
        acint.append( acint[it-1] + 2.*acov[it]/acov[0] )
    return acint


def tauintj(nt,nbins,ts,xhat = None):
    """ Given a time series, calculate jackknife bins of integrated autocorrelation time for each Markov time separation.

    INPUT:
          nt--The largest nt at which you think your estimate for tau_int could lie.
       nbins--The number of jackknife bins.
          ts--Time series array of measurements. Must be taken from equilibrium ensemble so that
              time translation invariance holds. List must be in order of markov chain generation
        xhat--True mean of time series (if you know it).

    OUTPUT:
      acintj--2D list indexed by time, then bin number acintj[it][ibin] """
    checkTS(ts)
    ndat=len(ts)
    if nbins<2:
        logger.TBError("Need nbins>1.")
    if nbins>=ndat:
        logger.TBError("Need nbins<ndat.")
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
        acor=[0.]*nbins
        # This is a check that we can't spill into the last bin
        itmax=nbins*binsize-binsize
        if it>=itmax:
            logger.TBError("it>=itmax.")
        for ibin in range(nbins):
            i1=ibin*binsize
            i2=(ibin+1)*binsize-1
            for i in range(i1,i2+1):
                acor[ibin]=acor[ibin]+(ts[i]-x)*(ts[i+it]-x)
            acor[ibin]=acor[ibin]/binsize
            # Bias correction
            if xhat is None:
                acor[ibin]=acor[ibin]*ndat/(ndat-1.)
        acorj.append(remove1Jackknife(acor))
    # Now make acintj
    acintj=[[1.]*nbins]
    for it in range(1,nt+1):
        tauintbins=[]
        for ibin in range(nbins):
            tauintbins.append(acintj[it-1][ibin]+2.*acorj[it][ibin]/acorj[0][ibin])
        acintj.append(tauintbins)
    return acintj


def getTauInt(ts, nbins, tpickMax, acoutfileName = 'acor.d', showPlot = False):
    """ Given a time series, return estimates for the integrated autocorrelation time and its error.

    INPUT:
         tpickMax--The largest nt where you think your estimate might become unreliable.
            nbins--The number of jackknife bins (for estimating the error in tau_int)
               ts--Time series array of measurements. Must be taken from equilibrium ensemble so that
                   time translation invariance holds. List must be in order of markov chain generation

    OUTPUT:
          tau_int--Estimate for integrated autocorrelation time.
         tau_inte--Its (jackknife) error bar.
           itpick--The Monte Carlo separation at which this method found its estimate for tau_int. """
    checkTS(ts)
    acoutfile=open(acoutfileName,'w')

    # Get integrated autocorrelation time list and corresponding jackknife list.
    acint  = np.array( tauint (tpickMax,ts,) )
    acintj = np.array( tauintj(tpickMax,nbins,ts) )

    # This block outputs time, tau_int and tau_int error, then gives an estimate of tau_int. When tau_int
    # decreases for the first time, we call that our estimate. This is because we know tau_int should be a 
    # monotonically increasing function of t.
    lmonoton=True
    tau_int=0.
    tau_inte=-1
    itpick=-1
    for it in range(tpickMax+1):
        ace = std_err(acintj[it])
        acoutfile.write(str(it)+'\t'+str(acint[it])+'\t'+str(ace)+'\n')
        if lmonoton:
            if not acint[it]<tau_int:
                tau_int=acint[it]
                tau_inte=ace
                itpick=it
            else:  # acint[it] < tau_int ==> tau_int decreased
                lmonoton=False
    acoutfile.close()

    if showPlot:
        clearPlot()
        plot_file(acoutfileName, xcol=0, ycol=1, yecol=2, xlabel='conf', ylabel='$\\tau_{\\rm int}$')
        plt.show()

    return tau_int, tau_inte, itpick
