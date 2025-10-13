# 
# lattice_params.py
# 
# D. Clarke
# 
# Class to handle input parameters of lattice configs. This is in particular for use with the HotQCD collaboration.
#
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import convert, fk_phys, r0_phys, r1_phys
from latqcdtools.physics.referenceScales import a_div_r1, a_times_fk, r0_div_a, CY_param, CY_phys
from latqcdtools.base.check import checkDomain
from latqcdtools.base.utilities import isReal


def _getMassString(mass):
    if mass is None:
        return None
    elif isinstance(mass,str):
        return mass
    elif isReal(mass):
        if not 0<mass<1:
            logger.TBRaise('I can only handle bare masses between 0 and 1.') 
        return str(mass)[2:] 
    else:
        logger.TBRaise(f'Expected str or float for mass. Got {type(mass)}')


def _getMassFloat(mass):
    if mass is None:
        return None
    elif isReal(mass):
        if not 0<mass<1:
            logger.TBRaise('I can only handle bare masses between 0 and 1.') 
        return mass
    elif isinstance(mass,str):
        fmass = float(f'0.{mass}')
        if not 0<fmass<1:
            logger.TBRaise('I can only handle bare masses between 0 and 1.') 
        return fmass
    else:
        logger.TBRaise(f'Expected str or float for mass. Got {type(mass)}')


class latticeParams:
    """
    A class to handle and check the input parameters of a lattice run.
    """

    def setCoupling(self,coupling):
        if isinstance(coupling, str):
            self.beta  = int(coupling)/10**(len(coupling)-1)
            self.cbeta = coupling
        else:
            self.beta  = coupling
            self.cbeta = str(coupling).replace('.','')

    def setScales(self,scaleYear,paramYear):
        self.fK    = fk_phys(CY_phys['fk'],"MeV")
        self.r1    = r1_phys(CY_phys['r1'],"fm")
        self.r0    = r0_phys(CY_phys['r0'],"fm")
        if paramYear is None:
            self.year = CY_param[self.scale] 
        else:
            self.year = paramYear
        if scaleYear is not None:
            if self.scale == 'fk':
                self.fK = fk_phys(scaleYear)
            elif self.scale == 'r1':
                self.r1 = r1_phys(scaleYear)
            elif self.scale == 'r0':
                self.r0 = r0_phys(scaleYear)

    #           mass1  mass2  mass3
    # Nf=1+1+1     mu     md     ms
    # Nf=2+1       ml     ms
    # Nf=X         ml   mpre
    # Nf=2+1+1     ml     ms     mc
    def __init__(self, Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None,
                 Nf='21', scaleYear=None, mu=0):
        """ 
        Based on some input, determine all parameters relevant to the ensemble.

        Parameters
        ----------
            Nsigma : int
            Ntau : int
            coupling : str, float
                beta = 6/g^2, where g is the bare coupling.
            mass1 : str
                Light-type mass.
            mass2 : str
                Strange-type mass. For degenerate Nf, can also be an auxiliary mass.
            mass3 : str
                Charm-type mass.
            scaleType : str
                Reference scale.
            paramYear : int
                Year during which the parameters of the rational function or polynomial function used for
                continuum limit extrapolations were determined.
            Nf : str
            scaleYear : int
                Year during which the value of the scale in physical units was determined.
            mu : float
                Baryon chemical potential.
        """
        checkDomain(scaleType,list(CY_phys.keys()))
        self.scale = scaleType
        self.setScales(scaleYear,paramYear)
        self.setCoupling(coupling)
        self.Nc   = 3
        self.mu   = mu
        self.Ns   = Nsigma
        self.Nt   = Ntau
        self.Nf   = Nf
        self.vol4 = self.Ns**3 * self.Nt
        self.cm1  = _getMassString(mass1)
        self.cm2  = _getMassString(mass2)
        self.cm3  = _getMassString(mass3)
        self.cml  = None 
        self.cmu  = None 
        self.cmd  = None 
        self.cms  = None 
        self.ml   = None 
        self.mu   = None 
        self.md   = None 
        self.ms   = None 
        self.cm   = None
        self.cpre = None
        self.m    = None
        self.pre  = None
        self.msml = None
        self.msmu = None
        self.mdmu = None
        if Nf == '21':
            if mass3 is not None:
                logger.TBRaise('Nf=2+1 expects only 2 mass parameters.')
            self.cml  = _getMassString(mass1)
            self.cms  = _getMassString(mass2)
            self.ml   = _getMassFloat(mass1)
            self.ms   = _getMassFloat(mass2)
        elif Nf == '111':
            self.cmu  = _getMassString(mass1) 
            self.cmd  = _getMassString(mass2) 
            self.cms  = _getMassString(mass3)
            self.mu   = _getMassFloat(mass1)
            self.md   = _getMassFloat(mass2)
            self.ms   = _getMassFloat(mass3)
        elif Nf == '211':
            self.cml  = _getMassString(mass1)
            self.cms  = _getMassString(mass2)
            self.cmc  = _getMassString(mass3)
            self.ml   = _getMassFloat(mass1)
            self.ms   = _getMassFloat(mass2)
            self.mc   = _getMassFloat(mass3)
        elif Nf=='3' or Nf=='5':
            if mass3 is not None:
                logger.TBRaise('Degenerate Nf expects only 2 mass parameters.')
            self.cm   = _getMassString(mass1)
            self.cpre = _getMassString(mass2)
            self.m    = _getMassFloat(mass1)
            self.pre  = _getMassFloat(mass2)
        else:
            logger.TBRaise("Unsupported Nf",Nf)
        if (self.ml is not None) and (self.ms is not None):
            self.msml=int(round(self.ms/self.ml))
        if (self.mu is not None) and (self.md is not None):
            self.mdmu=int(round(self.md/self.mu))
        if (self.mu is not None) and (self.ms is not None):
            self.msmu=int(round(self.ms/self.mu))
        if (self.beta<1.) or (10.<self.beta):
            logger.TBRaise("Invalid beta =",self.beta)


    def __repr__(self) -> str:
        return "latticeParams"


    def geta(self,units='fm'):
        if self.scale=='fk':
            return convert( a_times_fk(self.beta,self.year)/self.fK, 'MeVinv', units )
        elif self.scale=='r1':
            return convert( a_div_r1(self.beta,self.year)*self.r1, 'fm', units )
        elif self.scale=='r0':
            return convert( self.r0/r0_div_a(self.beta,self.year), 'fm', units )


    def getT(self,units='MeV') -> float:
        return convert( 1/(self.geta('fm')*self.Nt), 'fminv', units )


    def getLs(self,units='MeVinv') -> float:
        if self.Ns is not None:
            return convert( self.Ns*self.geta('fm'), 'fm', units )
        else:
            logger.TBRaise('Must specify Ns get get Ls.')


    # A nicely formatted summary of the lattice parameters.
    def paramSummary(self):
        logger.info()
        logger.info("Lattice parameter summary: ")
        if self.scale == 'fk':
            logger.info("    fK = ",round(self.fK*np.sqrt(2),2),"/sqrt(2) [MeV] ")
        elif self.scale == 'r1':
            logger.info("    r1 = ",self.r1,"[fm] ")
        elif self.scale == 'r0':
            logger.info("    r0 = ",round(self.r0,4),"[fm] ")
        if self.Ns is not None:
            logger.info("    Ns = ",self.Ns)
        logger.info("    Nt = ",self.Nt)
        if self.ml is not None:
            logger.info("    ml = ",self.ml)
        if self.mu is not None:
            logger.info("    mu = ",self.mu)
        if self.md is not None:
            logger.info("    md = ",self.md)
        if self.ms is not None:
            logger.info("    ms = ",self.ms)
        if self.m is not None:
            logger.info("     m = ",self.m)
        if self.pre is not None:
            logger.info("   pre = ",self.pre)
        if self.msml is not None: 
            logger.info(" ms/ml = ",self.msml)
        if self.msmu is not None: 
            logger.info(" ms/mu = ",self.msmu)
        if self.mdmu is not None: 
            logger.info(" md/mu = ",self.mdmu)
        logger.info("    T  = ",round(self.getT(),2), "[MeV]")
        logger.info("    a  = ",round(self.geta(),4), "[fm]")
        if self.Ns is not None:
            logger.info("    Ls = ",round(self.getLs(),4), "1/[MeV]")
        logger.info("  beta = ",self.beta)
        logger.info()
