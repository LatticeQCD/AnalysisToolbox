# 
# simulationManagement.py                                                               
# 
# D. Clarke
# 
# This is some Python code to help assist with managing lattice simulations. 
# 

import glob
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger
from latqcdtools.base.plotting import plt, plot_lines, plot_hist, set_params, clearPlot,\
    saveFigure
from latqcdtools.statistics.statistics import std_mean, std_err, checkTS, KSTest_1side
import scipy as sp


def countConfigurations(targetFolder,name,delimiter='.'):
    """
    Count the number of configurations in the target folder.

    Args:
        targetFolder (str)
        name (str): Assume configuration name has form name.###. 
        delimiter (str, optional): The delimiter between name and ###. Defaults to '.'.

    Returns:
        int: number of configurations in targetFolder 
    """
    checkType(str,targetFolder=targetFolder)
    checkType(str,name=name)
    checkType(str,delimiter=delimiter)
    confs=[]
    files=list(glob.iglob(f'{targetFolder}/{name}*'))
    for f in files:
        try:
            _, confno = f.rsplit(delimiter,1)
            confs.append(int(confno))
        except ValueError:
            logger.TBRaise(f'Configuration names in {targetFolder} must end in int. file={f}, delimiter="{delimiter}"')
    logger.info(f'{targetFolder}:  min={min(confs)}  max={max(confs)}')
    return len(confs)


def analyzeChain(MCtime,measurements,obslabel=None,MClabel=None,KScutoff=0.05,
                 showPlots=False,savePlots=False,plotNamePrefix=None,**plotargs):
    """
    Do some basic analysis of a MCMC time series of measurements. We check whether the data
    are distributed normally. Optionally you can plot the time series and/or a histogram.

    Args:
        MCtime (array-like): Array indexing the MC time. 
        measurements (array-like): Time series of measurements.
        obslabel (str, optional): Label for observable in plots. Defaults to None.
        MClabel (str, optional): Label for x-axis of time series plot. Defaults to None.
        KScutoff (float, optional): Threshold below which "not normal enough". Defaults to 0.05.
        showPlots (bool, optional): Defaults to False.
        savePlots (bool, optional): Defaults to False.
        plotNamePrefix (str, optional): Prefix of plot names. Defaults to None.
    """

    checkTS(MCtime)
    checkTS(measurements)
    checkType('real',KScutoff=KScutoff)

    if len(plotargs)>0 and showPlots==False and savePlots==False:
        logger.TBRaise('Passed plotting keywords without wanting to plot.')
    if savePlots and (plotNamePrefix is None):
        logger.TBRaise('Please specify a name prefix for the TS and histogram plots.')
    if savePlots or showPlots:
        if MClabel is None:
            logger.TBRaise('Plotting without an observable label.')
        checkType(str,MClabel=MClabel)

    plot_lines(MCtime,measurements,marker=None,color='black')
    set_params(xlabel=MClabel,ylabel=obslabel,**plotargs)
    if savePlots:
        checkType(str,plotNamePrefix=plotNamePrefix)
        saveFigure(plotNamePrefix+'TS.pdf')
    if showPlots:
        plt.show()
    clearPlot()

    plot_hist(measurements)
    set_params(xlabel=obslabel,**plotargs)
    if savePlots:
        saveFigure(plotNamePrefix+'hist.pdf')
    if showPlots:
        plt.show()
    clearPlot()

    mean = std_mean(measurements)
    err  = std_err(measurements)
    normalCDF = sp.stats.norm(loc=mean,scale=err).cdf

    if KSTest_1side(measurements,normalCDF)<KScutoff:
        logger.warn('Measurement distribution not consistent with Gaussian')
