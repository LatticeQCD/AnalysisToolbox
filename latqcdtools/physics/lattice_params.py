# 
# lattice_params.py
# 
# D. Clarke
# 
# Class to handle input parameters of lattice configs. This is in particular for use with the HotQCD collaboration.
#
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import MeVinv_to_fm, fm_to_MeVinv, fk_phys, r0_phys, r1_phys
from latqcdtools.physics.referenceScales import a_div_r1, a_times_fk, r0_div_a, CY_param, CY_phys
from latqcdtools.base.check import checkDomain


def _massStringToFloat(string):
    if string is None:
        return None
    else:
        return float('0.' + string)


class latticeParams:
    """A class to handle and check the input parameters of a lattice run."""

    # If doing Nf=2+1 physics, we interpret mass1 and mass2 as light and strange masses, respectively. If doing
    # degenerate Nf physics, we interpret mass1 and mass2 as the quark mass and preconditioner, respectively.
    def __init__(self, Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None,
                 Nf='21', scaleYear=None, mu=0):
        """ Based on some input, determine all parameters relevant to the ensemble.

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
        checkDomain(scaleType,CY_phys.keys())
        self.scale = scaleType
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
        if isinstance(coupling, str):
            self.beta  = int(coupling)/10**(len(coupling)-1)
            self.cbeta = coupling
        else:
            self.beta  = coupling
            self.cbeta = str(coupling).replace('.','')
        self.Nc   = 3
        self.mu   = mu
        self.cm1  = mass1
        self.cm2  = mass2
        self.cm3  = mass3
        self.Ns   = Nsigma
        self.Nt   = Ntau
        self.Nf   = Nf
        self.vol4 = self.Ns**3 * self.Nt
        if Nf=='21':
            if mass3 is not None:
                logger.TBError('Nf=2+1 expects only 2 mass parameters.')
            self.cml  = mass1
            self.cms  = mass2
            self.ml   = _massStringToFloat(mass1)
            self.ms   = _massStringToFloat(mass2)
            self.cm   = None
            self.cpre = None
            self.m    = None
            self.pre  = None
        elif Nf == '211':
            self.cml  = mass1
            self.cms  = mass2
            self.cmc  = mass3
            self.ml   = _massStringToFloat(mass1)
            self.ms   = _massStringToFloat(mass2)
            self.mc   = _massStringToFloat(mass3)
            self.cm   = None
            self.cpre = None
            self.m    = None
            self.pre  = None
        elif Nf=='3' or Nf=='5':
            if mass3 is not None:
                logger.TBError('Degenerate Nf expects only 2 mass parameters.')
            self.cml  = None
            self.cms  = None
            self.ml   = None
            self.ms   = None
            self.cm   = mass1
            self.cpre = mass2
            self.m    = _massStringToFloat(mass1)
            self.pre  = _massStringToFloat(mass2)
        if (self.ml is not None) and (self.ms is not None):
            self.msml=int(round(self.ms/self.ml))
        if (self.beta<1.) or (10.<self.beta):
            logger.TBError("Invalid beta =",self.beta)
        if (scaleType!='fk') and (scaleType!='r0') and (scaleType!='r1'):
            logger.TBError("Unknown reference scale",scaleType)


    def __repr__(self) -> str:
        return "latticeParams"


    # a in [fm]
    def geta(self):
        if self.scale=='fk':
            return MeVinv_to_fm( a_times_fk(self.beta,self.year)/self.fK )
        elif self.scale=='r1':
            return a_div_r1(self.beta,self.year)*self.r1
        elif self.scale=='r0':
            return self.r0/r0_div_a(self.beta,self.year)


    def getT(self) -> float:
        """ T in MeV. """
        return 1/fm_to_MeVinv( (self.geta()*self.Nt) )


    def getLs(self) -> float:
        """ L in space-like direction in [1/MeV]. """
        if self.Ns is not None:
            return fm_to_MeVinv(self.Ns*self.geta())
        else:
            logger.TBError('Must specify Ns get get Ls.')


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
        if self.ms is not None:
            logger.info("    ms = ",self.ms)
        if self.m is not None:
            logger.info("     m = ",self.m)
        if self.pre is not None:
            logger.info("   pre = ",self.pre)
        if (self.ml is not None) and (self.ms is not None): 
            logger.info(" ms/ml = ",self.msml)
        logger.info("    T  = ",round(self.getT(),2), "[MeV]")
        logger.info("    a  = ",round(self.geta(),4), "[fm]")
        if self.Ns is not None:
            logger.info("    Ls = ",round(self.getLs(),4), "1/[MeV]")
        logger.info("  beta = ",self.beta)
        logger.info()
