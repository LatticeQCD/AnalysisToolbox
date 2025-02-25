latqcdtools.interfaces.simulationManagement
=============

`analyzeChain(MCtime, measurements, obslabel=None, MClabel=None, KScutoff=0.05, showPlots=False, savePlots=False, plotNamePrefix=None, **plotargs)`

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

`countConfigurations(targetFolder, name, delimiter='.')`

Count the number of configurations in the target folder.

Args:
    targetFolder (str)
    name (str): Assume configuration name has form name.###. 
    delimiter (str, optional): The delimiter between name and ###. Defaults to '.'.

Returns:
    int: number of configurations in targetFolder 

