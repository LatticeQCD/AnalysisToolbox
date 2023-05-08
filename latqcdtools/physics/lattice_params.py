# 
# lattice_params.py
# 
# D. Clarke
# 
# Class to handle input parameters of lattice configs. This is in particular for use with the HotQCD collaboration.
#
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import MeVinv_to_fm, fm_to_MeVinv
from latqcdtools.physics.referenceScales import fk_phys, r1_MILC_2010, a_div_r1, a_times_fk, r0_div_a, r0_hQCD_2014, CY_A_TIMES_FK, CY_FK_PHYS


def massStringToFloat(string):
    if string is None:
        return None
    else:
        return float('0.' + string)


class latticeParams:
    """A class to handle and check the input parameters of a lattice run."""

    r1=r1_MILC_2010("fm")
    r0=r0_hQCD_2014("fm")

    # If doing Nf=2+1 physics, we interpret mass1 and mass2 as light and strange masses, respectively. If doing
    # degenerate Nf physics, we interpret mass1 and mass2 as the quark mass and preconditioner, respectively.
    def __init__(self, Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=CY_A_TIMES_FK,
                 Nf='21', scaleYear=CY_FK_PHYS, mu=0):
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
        self.fK = fk_phys(scaleYear,"MeV")
        if isinstance(coupling, str):
            self.beta  = int(coupling)/10**(len(coupling)-1)
            self.cbeta = coupling
        else:
            self.beta  = coupling
            self.cbeta = str(coupling).replace('.','')
        self.Nc    = 3
        self.mu    = mu
        self.cm1   = mass1
        self.cm2   = mass2
        self.cm3   = mass3
        self.year  = paramYear
        self.Ns    = Nsigma
        self.Nt    = Ntau
        self.scale = scaleType
        self.Nf    = Nf
        self.vol4  = self.Ns**3 * self.Nt
        if Nf=='21':
            if mass3 is not None:
                logger.TBError('Nf=2+1 expects only 2 mass parameters.')
            self.cml  = mass1
            self.cms  = mass2
            self.ml   = massStringToFloat(mass1)
            self.ms   = massStringToFloat(mass2)
            self.cm   = None
            self.cpre = None
            self.m    = None
            self.pre  = None
        elif Nf == '211':
            self.cml  = mass1
            self.cms  = mass2
            self.cmc  = mass3
            self.ml   = massStringToFloat(mass1)
            self.ms   = massStringToFloat(mass2)
            self.mc   = massStringToFloat(mass3)
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
            self.m    = massStringToFloat(mass1)
            self.pre  = massStringToFloat(mass2)
        if (self.ml is not None) and (self.ms is not None):
            self.msml=int(round(self.ms/self.ml))
        if (self.beta<1.) or (10.<self.beta):
            logger.TBError("Invalid beta =",self.beta)
        if (scaleType!='fk') and (scaleType!='r0') and (scaleType!='r1'):
            logger.TBError("Unknown reference scale",scaleType)


    # a in [fm]
    def geta(self):
        if self.scale=='fk':
            return MeVinv_to_fm( a_times_fk(self.beta,self.year)/self.fK )
        elif self.scale=='r1':
            return a_div_r1(self.beta,self.year)*self.r1
        elif self.scale=='r0':
            return self.r0/r0_div_a(self.beta,self.year)


    # T in [MeV]
    def getT(self):
        if self.Ns is not None:
            if self.Ns == self.Nt:
                return 0.
        return 1/fm_to_MeVinv( (self.geta()*self.Nt) )
    
    
    # L in spacelike direction in [1/MeV]
    def getLs(self):
        if self.Ns is not None:
            return fm_to_MeVinv(self.Ns*self.geta())
        else:
            logger.TBError('Must specify Ns get get Ls.')


    # A nicely formatted summary of the lattice parameters.
    def paramSummary(self):
        print("\nLattice parameter summary: ")
        if self.scale == 'fk':
            print("    fK = ",round(self.fK*np.sqrt(2),2),"/sqrt(2) [MeV] ")
        elif self.scale == 'r1':
            print("    r1 = ",self.r1,"[fm] ")
        elif self.scale == 'r0':
            print("    r0 = ",round(self.r0,4),"[fm] ")
        if self.Ns is not None:
            print("    Ns = ",self.Ns)
        print("    Nt = ",self.Nt)
        if self.ml is not None:
            print("    ml = ",self.ml)
        if self.ms is not None:
            print("    ms = ",self.ms)
        if self.m is not None:
            print("     m = ",self.m)
        if self.pre is not None:
            print("   pre = ",self.pre)
        if (self.ml is not None) and (self.ms is not None): 
            print(" ms/ml = ",self.msml)
        print("    T  = ",round(self.getT(),2), "[MeV]")
        print("    a  = ",round(self.geta(),4), "[fm]")
        if self.Ns is not None:
            print("    Ls = ",round(self.getLs(),4), "1/[MeV]")
        print("  beta = ",self.beta,"\n")
