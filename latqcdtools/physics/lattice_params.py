# 
# lattice_params.py
# 
# D. Clarke
# 
# Class to handle input parameters of lattice configs. This is in particular for use with the HotQCD collaboration.
#
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.unitConversions import MeVinv_to_fm, fm_to_MeVinv
from latqcdtools.physics.referenceScales import fk_PDG_2012, r1_MILC_2010, a_div_r1, a_times_fk, r0_div_a, r0_hQCD_2014


def massStringToFloat(string):
    if string is None:
        return None
    else:
        return float('0.' + string)


class latticeParams:
    """A class to handle and check the input parameters of a lattice run."""

    fK=fk_PDG_2012("MeV")
    r1=r1_MILC_2010("fm")
    r0=r0_hQCD_2014("fm")

    # If doing Nf=2+1 physics, we interpret mass1 and mass2 as light and strange masses, respectively. If doing
    # degenerate Nf physics, we interpret mass1 and mass2 as the quark mass and preconditioner, respectively.
    def __init__(self, Nsigma, Ntau, coupling, mass1=None, mass2=None, scaleType='fk', paramYear=2021, Nf='21'):
        if isinstance(coupling, str):
            self.beta  = int(coupling)/1000
            self.cbeta = coupling
        else:
            self.beta  = coupling
            self.cbeta = str(int(coupling*1000))
        self.cm1   = mass1
        self.cm2   = mass2
        self.year  = paramYear
        self.Ns    = Nsigma
        self.Nt    = Ntau
        self.vol4  = self.Ns**3 * self.Nt
        self.vol3  = self.Ns**3
        self.scale = scaleType
        self.Nf    = Nf
        if Nf=='21':
            self.cml  = mass1
            self.cms  = mass2
            self.ml   = massStringToFloat(mass1)
            self.ms   = massStringToFloat(mass2)
            self.cm   = None
            self.cpre = None
            self.m    = None
            self.pre  = None
        else:
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
            logger.TBError("Invalid beta.")
        if (scaleType!='fk') and (scaleType!='r0') and (scaleType!='r1'):
            logger.TBError("Unknown reference scale",scaleType)


    # a in [fm]
    def geta(self):
        if self.scale=='fk':
            return MeVinv_to_fm( a_times_fk(self.beta,self.year)/self.fK )
        elif self.scale=='r1':
            return a_div_r1(self.beta,self.year)*self.r1
        elif self.scale=='r0':
            return self.r0/r0_div_a(self.beta)


    # T in [MeV]
    def getT(self):
        if self.Ns == self.Nt:
            return 0.
        else:
            return 1/fm_to_MeVinv( (self.geta()*self.Nt) )


    # A nicely formatted summary of the lattice parameters.
    def paramSummary(self):
        print("\nLattice parameter summary: ")
        if self.scale == 'fk':
            print("    fK = ",round(self.fK*np.sqrt(2),2),"/sqrt(2) [MeV] ")
        elif self.scale == 'r1':
            print("    r1 = ",self.r1,"[fm] ")
        elif self.scale == 'r0':
            print("    r0 = ",round(self.r0,4),"[fm] ")
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
        print("  beta = ",self.beta,"\n")
