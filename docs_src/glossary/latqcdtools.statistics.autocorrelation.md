latqcdtools.statistics.autocorrelation
=============

```python
getTauInt(ts, nbins, tpickMax, acoutfileName=None, showPlot=False):
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
```
```python
gewekeZ(ts, first=0.1, last=0.5, nbins=5) -> float:
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
```
```python
remove1Jackknife(ts) -> numpy.ndarray:
"""
Create remove-1 jackknife list from 1-d series.

Args:
    ts (array-like): time series 

Returns:
    np.array: 1-d array of jackknife means 
"""
```
```python
tauint(nt, ts, xhat=None) -> numpy.ndarray:
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
```
```python
tauintj(nt, nbins, ts, xhat=None) -> numpy.ndarray:
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
```
