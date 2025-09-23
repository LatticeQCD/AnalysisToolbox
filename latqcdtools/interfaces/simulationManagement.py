# 
# simulationManagement.py                                                               
# 
# D. Clarke
# 
# This is some Python code to help assist with managing lattice simulations. 
# 

from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger
from latqcdtools.base.plotting import plt, plot_lines, plot_hist, set_params, clearPlot,\
    saveFigure
from latqcdtools.statistics.statistics import std_mean, std_err, checkTS, KSTest_1side
import scipy as sp
from latqcdtools.base.utilities import ls, naturalSort
from latqcdtools.math.math import rel_check

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
    files=ls(f'{targetFolder}/{name}*')
    for f in files:
        try:
            _, confno = f.rsplit(delimiter,1)
            confs.append(int(confno))
        except ValueError as e:
            logger.TBRaise(f'Configuration names in {targetFolder} must end in int. file={f}, delimiter="{delimiter}"',exception=e)
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

    if savePlots or showPlots:
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



_allowedEnsKinds = ['conf','meas','conf/meas','seed','seed/meas']
_allowedEnsKeys  = ['kind',  # configurations? measurements?
                    'date',  # last date you checked it
                    'Nt',
                    'Ns', 
                    'beta',
                    'Nf', 
                    'ml',    # Nf=2+1: light mass 
                    'ms',    # Nf=2+1: strange mass
                    'mc',    # Nf=2+1: charm mass
                    'mf',    # Nf=N: fermion mass
                    'mpre',  # Nf=N: preconditioner mass
                    'location', 
                    'folder', 
                    'notes', 
                    'creator'
                    ]
_allowedEnsNf    = ['2+1','3','5','2+1+1']


class ensemble:

    """ 
    An ensemble object. This is to help you organize where data are saved
    """

    def __repr__(self) -> str:
        return "ensemble"

    def __init__(self,**kwargs):
        self._ens = {
            'kind'    :None,
            'date'    :None, 
            'Nt'      :None, 
            'Ns'      :None, 
            'beta'    :None, 
            'Nf'      :None,
            'ml'      :None, 
            'ms'      :None, 
            'mc'      :None, 
            'mf'      :None, 
            'mpre'    :None,
            'ml/ms'   :None, 
            'location':None, 
            'folder'  :None, 
            'notes'   :None,
            'creator' :None
        }
        for key in kwargs:
            if key not in _allowedEnsKeys:
                logger.warn(f'Unknown ensemble keyword {key}')
            else:
                self._ens[key] = kwargs[key] 
        lcomplete = True
        if self._ens['kind'] is None:
            lcomplete = False
        if self._ens['Nt'] is None:
            lcomplete = False
        if self._ens['Nf'] is None:
            lcomplete = False
        if self._ens['beta'] is None:
            lcomplete = False
        if not lcomplete:
            logger.TBRaise('Must minimally specify Nt, Nf, beta, and kind.')
        if self._ens['kind'] not in _allowedEnsKinds: 
            logger.TBRaise('Unrecognized kind',self._ens['kind'])
        if self._ens['Nf'] not in _allowedEnsNf: 
            logger.TBRaise('Unrecognized Nf',self._ens['Nf'])
        if self._ens['Nf']=='2+1' or self._ens['Nf']=='2+1+1':
            if (self._ens['mf'] is not None) or (self._ens['mpre'] is not None):
                logger.TBRaise('Detected degenerate Nf param with Nf=2+1(+1)')
        if self._ens['Nf']=='3' or self._ens['Nf']=='5':
            if (self._ens['ms'] is not None) or (self._ens['ms'] is not None) or (self._ens['mc'] is not None):
                logger.TBRaise('Detected Nf=2+1(+1) param with degenerate Nf')
        if self._ens['mc'] is not None:
            if self._ens['Nf']!='2+1+1':
                logger.TBRaise('Specified charm mass without Nf=2+1+1')
        if (self._ens['ms'] is not None) and (self._ens['ml'] is not None):
            msml = self._ens['ms']/self._ens['ml']
            if rel_check(msml,27,0.01):
                self._ens['ml/ms'] = 'phys'
            if rel_check(msml,80,0.01):
                self._ens['ml/ms'] = '1/80'

    def get(self,key):
        if self._ens[key] is not None:
            return self._ens[key]
        else:
            return '-'

    def accessor(self) -> dict:
        return self._ens


class repository():

    """
    This reposistory holds ensemble objects.
    """

    def __init__(self):
        self._repo = []

    def accessor(self) -> list:
        return self._repo

    def append(self,ens):
        checkType(ensemble,ens=ens) 
        self._repo.append(ens)

    def search(self,key,val,verbose=True):
        matches = repository()
        nmatches = 0
        for ens in self._repo: 
            if ens.get(key)==val:
                matches.append(ens)
                nmatches += 1
        if verbose:
            logger.info(f'Found {nmatches} matches:')
            matches.list(verbose=True)
        return matches

    def list(self,verbose=False):
        if verbose:
            logger.info()
            for ens in self._repo:
                for key in _allowedEnsKeys: 
                    if ens.get(key) is not None:
                        logger.info(f'  {key} = {ens.get(key)}')
                logger.info()
        logger.info(f'N_ens = {len(self._repo)}')

    def findUniqueKeys(self,key) -> set:
        ret = []
        for ens in self._repo:
            ret.append(ens.get(key))
        return set(ret)

    def table(self):
        logger.info(f'{"Nf":<6} {"ml/ms":<6} {"Nt":<4} {"beta":<8} {"Ns":<4} {"ml (mf)":<12} {"ms (mpre)":<10} {"mc":<8} {"kind":<12} {"location":<10} folder')
        logger.info()

        # This loop controls the order in which entries appear in the table.
        Nfs = self.findUniqueKeys('Nf') 
        for Nf in naturalSort(Nfs):
            repoNf = self.search('Nf',Nf,verbose=False)
            mqs = repoNf.findUniqueKeys('ml/ms')
            for mq in mqs:
                repomq = repoNf.search('ml/ms',mq,verbose=False)
                Nts = repomq.findUniqueKeys('Nt') 
                for Nt in sorted(Nts):
                    repoNt = repomq.search('Nt',Nt,verbose=False)
                    betas = repoNt.findUniqueKeys('beta')
                    for beta in sorted(betas):
                        repoBeta = repoNt.search('beta',beta,verbose=False)
                        Nss = repoBeta.findUniqueKeys('Ns')
                        for Ns in sorted(Nss):
                            repoNs = repoBeta.search('Ns',Ns,verbose=False)
                            for ens in repoNs.accessor():
                                if Nf=='2+1':
                                    logger.info(f"{Nf:<6} {ens.get('ml/ms'):<6} {Nt:<4} {beta:<8} {Ns:<4} {ens.get('ml'):<12} {ens.get('ms'):<10} {'-':<8} {ens.get('kind'):<12} {ens.get('location'):<10} {ens.get('folder')}")
                                elif Nf=='2+1+1':
                                    logger.info(f"{Nf:<6} {ens.get('ml/ms'):<6} {Nt:<4} {beta:<8} {Ns:<4} {ens.get('ml'):<12} {ens.get('ms'):<10} {ens.get('mc'):<8} {ens.get('kind'):<12} {ens.get('location'):<10} {ens.get('folder')}")
                                elif Nf=='3' or Nf=='5':
                                    logger.info(f"{Nf:<6} {ens.get('ml/ms'):<6} {Nt:<4} {beta:<8} {Ns:<4} {ens.get('mf'):<12} {ens.get('mpre'):<10} {'-':<8} {ens.get('kind'):<12} {ens.get('location'):<10} {ens.get('folder')}")