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
from latqcdtools.physics.scales_hisq import fk_PDG_2012, r1_MILC_2010, a_div_r1, a_times_fk
from latqcdtools.physics.scales_quenched import r0_div_a, r0_hQCD_2014


def massRatioToMasses(msml, Nt, cbeta):
    if Nt==8 and msml==160:
        massTable=quarkMassNt8.Table160
    elif Nt==8 and msml==80:
        massTable=quarkMassNt8.Table80
    elif Nt==8 and msml==40:
        massTable=quarkMassNt8.Table40
    elif Nt==8 and msml==27:
        massTable=quarkMassNt8.Table27
    elif Nt==12 and msml==27:
        massTable=quarkMassNt12.Table27
    elif msml is None:
        pass
    else:
        logger.TBError("ms/ml not correctly set.")
    if msml is None:
        cml=None
        cms=None
    else:
        cml=massTable[cbeta][0]
        cms=massTable[cbeta][1]
    return cml, cms


def massStringToFloat(string):
    if string is None:
        return -1
    else:
        return float('0.' + string)


class latticeParams:
    """A class to handle and check the input parameters of a lattice run."""

    fK=fk_PDG_2012("MeV")
    r1=r1_MILC_2010("fm")
    r0=r0_hQCD_2014("fm")

    def __init__(self, Nsigma, Ntau, coupling, mass_l=None, mass_s=None, scaleType='fk', paramYear=2021):
        if isinstance(coupling, str):
            self.beta  = int(coupling)/1000
            self.cbeta = coupling
        else:
            self.beta  = coupling
            self.cbeta = str(int(coupling*1000))
        self.year  = paramYear
        self.Ns    = Nsigma
        self.Nt    = Ntau
        self.cml   = mass_l
        self.cms   = mass_s
        self.ml    = massStringToFloat(mass_l)
        self.ms    = massStringToFloat(mass_s)
        self.vol4  = self.Ns**3 * self.Nt
        self.vol3  = self.Ns**3
        self.scale = scaleType
        if (self.ml is not None) and (self.ms is not None):
            self.msml=int(round(self.ms/self.ml))
        if (self.scale=='r0') and ( (mass_l is not None) or (mass_s is not None) ):
            logger.warn("Using pure SU(3) scale for 2+1 flavor QCD.")
        if ( (self.scale=='r1') or (self.scale=='fk') ) and (mass_l is None) and (mass_s is None) :
            logger.warn("Using 2+1 flavor QCD scale for pure SU(3).")
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
        if self.ml >= 0.:
            print("    ml = ",self.ml)
        if self.ms >= 0.:
            print("    ms = ",self.ms)
        if (self.ml > 0.) and (self.ms >= 0.):
            print(" ms/ml = ",self.msml)
        print("    T  = ",round(self.getT(),2), "[MeV]")
        print("    a  = ",round(self.geta(),4), "[fm]")
        print("  beta = ",self.beta,"\n")


    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        return self.getcgeom()+'f21b'+self.cbeta+'m'+self.cml+'m'+self.cms
    def getcGradFlowPureGauge(self):
        return 's'+str(self.Ns).zfill(3)+'t'+str(self.Nt).zfill(2)+'_b0'+self.cbeta+'00'


class quarkMassNt12:

    """Lookup tables for Nt=12, Nf=2+1 flavors. HotQCD HISQ runs."""

    Table27={ '6794': ['00167','0450'],
              '6850': ['00157','0424'],
              '6910': ['00148','0401'] }


class quarkMassNt8:

    """Lookup tables for Nt=8, Nf=2+1 flavors. HotQCD HISQ runs."""

    Table160={ '6285': ['00049375','0790'],
               '6300': ['00048250','0772'],
               '6315': ['00047500','0760'],
               '6330': ['00046625','0746'],
               '6354': ['00045500','0728'],
               '6372': ['00044456','0711'],
               '6390': ['00043375','0694'],
               '6423': ['00041875','0670'],
               '6445': ['00040750','0652'] }

    Table80={ '6285': ['0009875','0790'],
              '6300': ['0009650','0772'],
              '6315': ['0009500','0760'],
              '6330': ['0009325','0746'],
              '6354': ['0009100','0728'],
              '6372': ['0008891','0711'],
              '6390': ['0008675','0694'],
              '6423': ['0008375','0670'],
              '6445': ['0008150','0652'] }

    Table40={ '6260': ['002025','0810'],
              '6285': ['001975','0790'],
              '6300': ['001930','0772'],
              '6315': ['001900','0760'],
              '6330': ['001865','0746'],
              '6354': ['001820','0728'],
              '6365': ['001790','0716'],
              '6390': ['001735','0694'],
              '6423': ['001675','0670'],
              '6445': ['001630','0652'],
              '6474': ['001580','0632'],
              '6500': ['001535','0614'] }

    Table27={ '6245': ['00307','0830'],
              '6285': ['00293','0790'],
              '6315': ['00281','0759'],
              '6354': ['00270','0728'],
              '6390': ['00257','0694'],
              '6423': ['00248','0670'],
              '6445': ['00241','0652'],
              '6474': ['00234','0632'],
              '6500': ['00228','0614'] }
