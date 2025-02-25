latqcdtools.statistics.autocorrelation
=============

`getTauInt(ts, nbins, tpickMax, acoutfileName='acor.d', showPlot=False)`

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

`remove1Jackknife(ts) -> numpy.ndarray`

Create remove-1 jackknife list from 1-d series.

Args:
    ts (array-like): time series 

Returns:
    np.array: 1-d array of jackknife means 

`tauint(nt, ts, xhat=None) -> numpy.ndarray`

Given a time series, calculate estimators for its integrated autocorrelation time  at each Markov time separation.

INPUT:
     nt--The largest you think tau_int could be.
     ts--Time series array of measurments. Must be taken from equilibrium ensemble so that
         time translation invariance holds. List must be in order of Markov chain generation.
   xhat--True mean of time series (if you know it).

OUTPUT:
  acint--List of integrated autocorrelation times. 
  
`tauintj(nt, nbins, ts, xhat=None) -> numpy.ndarray`

Given a time series, calculate jackknife bins of integrated autocorrelation time for each Markov time separation.

INPUT:
      nt--The largest nt at which you think your estimate for tau_int could lie.
   nbins--The number of jackknife bins.
      ts--Time series array of measurements. Must be taken from equilibrium ensemble so that
          time translation invariance holds. List must be in order of markov chain generation
    xhat--True mean of time series (if you know it).

OUTPUT:
  acintj--2D list indexed by time, then bin number acintj[it][ibin] 
  
