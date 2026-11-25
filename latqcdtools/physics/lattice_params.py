# 
# lattice_params.py
# 
# D. Clarke, K. Ebira, Google Antigravity
# 
# Class to handle input parameters of lattice configs. This is in particular for use with the HotQCD collaboration.
#
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import convert, fk_phys, r0_phys, r1_phys
from latqcdtools.physics.referenceScales import a_div_r1, a_times_fk, r0_div_a, CY_param, CY_phys, get_aWorld
from latqcdtools.base.check import checkDomain
from latqcdtools.base.utilities import isReal, isIntType


CHECKWORLD = True
_FK_PHYS_WORLDS = ['nature', 'Nf2', 'Nf21', 'Nf211']
_FK_NATURE_YEARS = (2012, 2018)


def _parseYear(year, name):
    """
    Validate that a year argument is an integer (or integer string) and return it as int.
    """
    if isIntType(year):
        return int(year)
    if isinstance(year, str):
        stripped = year.strip()
        if stripped.isdigit():
            return int(stripped)
    logger.TBRaise(f"Expected int for {name}, got {type(year).__name__}: {year!r}")


def ignoreWorldWarning():
    """ 
    Turn off world warnings. 
    """
    global CHECKWORLD 
    CHECKWORLD = False
    logger.warn("Ignoring world warnings.")


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

    def setScales(self, scaleYear, paramYear, physWorld=None):
        """
        Determine the a(beta) parameterization world (aWorld), physical-units world (physWorld),
        and physical reference scale value for this ensemble.
        """
        self.fK = None
        self.r1 = None
        self.r0 = None
        self.physWorld = None
        self.scaleYear = None

        ensemble_world = 'Nf0' if self.Nf is None else f'Nf{self.Nf}'
        # Keep self.world backwards-compatible with constants.worlds ('SU3' for pure gauge) while using 'Nf0' for aWorld.
        self.world = 'SU3' if self.Nf is None else f'Nf{self.Nf}'

        # Smart paramYear default: Nf=2+1 r0 uses the 2012 2+1-flavor fit; pure gauge uses CY_param['r0'] (2017).
        if paramYear is None:
            if self.scale == 'r0':
                if self.Nf == '21':
                    self.year = 2012
                else:
                    self.year = CY_param['r0']
            else:
                self.year = CY_param[self.scale]
        else:
            self.year = _parseYear(paramYear, 'paramYear')

        self.aWorld = get_aWorld(self.scale, self.year)

        if CHECKWORLD and (ensemble_world != self.aWorld):
            logger.warn(f"Ensemble world '{ensemble_world}' does not match parameterization world '{self.aWorld}'.")

        parsed_scale_year = None if scaleYear is None else _parseYear(scaleYear, 'scaleYear')

        if self.scale == 'fk':
            if physWorld is None:
                if parsed_scale_year in _FK_NATURE_YEARS:
                    self.physWorld = 'nature'
                else:
                    self.physWorld = 'Nf21'
            else:
                checkDomain(physWorld, _FK_PHYS_WORLDS)
                self.physWorld = physWorld
            if parsed_scale_year is None:
                year = 2018 if self.physWorld == 'nature' else CY_phys['fk']
            else:
                year = parsed_scale_year
            self.scaleYear = year
            self.fK = fk_phys(year=year, units="MeV", returnErr=False, world=self.physWorld)
        elif self.scale == 'r1':
            if physWorld is not None and physWorld != 'Nf21':
                logger.TBRaise(f"Unsupported physWorld '{physWorld}' for scale 'r1'. Only 'Nf21' is available.")
            self.physWorld = 'Nf21'
            year = CY_phys['r1'] if parsed_scale_year is None else parsed_scale_year
            self.scaleYear = year
            self.r1 = r1_phys(year=year, units="fm", returnErr=False)
        elif self.scale == 'r0':
            if physWorld is not None and physWorld != 'Nf21':
                logger.TBRaise(f"Unsupported physWorld '{physWorld}' for scale 'r0'. Only 'Nf21' is available.")
            self.physWorld = 'Nf21'
            year = CY_phys['r0'] if parsed_scale_year is None else parsed_scale_year
            self.scaleYear = year
            self.r0 = r0_phys(year=year, units="fm", returnErr=False)
        else:
            if physWorld is not None:
                logger.TBRaise(f"physWorld is not yet supported for scale '{self.scale}'.")
            if parsed_scale_year is not None and parsed_scale_year != CY_phys[self.scale]:
                logger.TBRaise(f"Unsupported scaleYear {parsed_scale_year} for scale '{self.scale}'.")
            self.physWorld = None
            self.scaleYear = None

    #           mass1  mass2  mass3
    # Nf=1+1+1     mu     md     ms
    # Nf=2+1       ml     ms
    # Nf=X         ml   mpre
    # Nf=2+1+1     ml     ms     mc
    def __init__(self, Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None,
                 Nf=None, scaleYear=None, muB=0, physWorld=None):
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
                Reference scale ('fk', 'r1', or 'r0'; 't0' is accepted for metadata only, without geta/getT).
            paramYear : int
                Year during which the parameters of the rational function or polynomial function used for
                continuum limit extrapolations were determined.
            Nf : str or None
                Sea quark content ('21', '111', '211', '2', '3', '5', or None for pure gauge).
            scaleYear : int
                Year during which the value of the scale in physical units was determined.
            muB : float
                Baryon chemical potential.
            physWorld : str, optional
                World from which the physical-units scale value is taken. Selectable for scaleType='fk'
                ('nature', 'Nf2', 'Nf21', or 'Nf211'). When None, defaults to 'nature' for PDG scaleYear
                2012/2018 and 'Nf21' otherwise; for 'r0' and 'r1' only 'Nf21' is available.
        """
        checkDomain(scaleType,list(CY_phys.keys()))
        self.scale = scaleType
        self.Nf    = Nf
        self.Ns   = Nsigma
        self.Nt   = Ntau
        self.vol4 = (self.Ns**3 * self.Nt) if (self.Ns is not None and self.Nt is not None) else None
        self.setCoupling(coupling)
        self.Nc   = 3
        self.muB  = muB
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
        elif Nf in ('2', '3', '5'):
            if mass3 is not None:
                logger.TBRaise('Degenerate Nf expects only 2 mass parameters.')
            self.cm   = _getMassString(mass1)
            self.cpre = _getMassString(mass2)
            self.m    = _getMassFloat(mass1)
            self.pre  = _getMassFloat(mass2)
        elif Nf is None:
            if (mass1 is not None) or (mass2 is not None) or (mass3 is not None):
                logger.TBRaise('Nf=None (pure gauge) expects no mass parameters.')
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
        # Set scales after Nf and beta validation so invalid inputs fail cleanly without spurious world warnings.
        self.setScales(scaleYear, paramYear, physWorld=physWorld)


    def __repr__(self) -> str:
        Ns = getattr(self, 'Ns', None)
        Nt = getattr(self, 'Nt', None)
        beta = getattr(self, 'beta', None)
        scale = getattr(self, 'scale', None)
        return f"latticeParams(Ns={Ns}, Nt={Nt}, beta={beta}, scale={scale})"


    def geta(self,units='fm'):
        if self.scale=='fk':
            return convert( a_times_fk(self.beta,self.year)/self.fK, 'MeVinv', units )
        elif self.scale=='r1':
            return convert( a_div_r1(self.beta,self.year)*self.r1, 'fm', units )
        elif self.scale=='r0':
            return convert( self.r0/r0_div_a(self.beta,self.year), 'fm', units )
        else:
            logger.TBRaise(f"Lattice spacing extraction not yet implemented for scale '{self.scale}'.")


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
        logger.info("       Nf =", self.Nf)
        logger.info("   aWorld =", f"{self.aWorld} ({self.scale} {self.year})")
        if self.physWorld is not None:
            logger.info("physWorld =", f"{self.physWorld} ({self.scaleYear})")
        else:
            logger.info("physWorld =", self.physWorld)
        if self.scale == 'fk':
            logger.info("       fK =", round(self.fK*np.sqrt(2), 2), "/sqrt(2) [MeV] ")
        elif self.scale == 'r1':
            logger.info("       r1 =", self.r1, "[fm] ")
        elif self.scale == 'r0':
            logger.info("       r0 =", round(self.r0, 4), "[fm] ")
        if self.Ns is not None:
            logger.info("       Ns =", self.Ns)
        if self.Nt is not None:
            logger.info("       Nt =", self.Nt)
        logger.info("     beta =", self.beta)
        if self.ml is not None:
            logger.info("       ml =", self.ml)
        if self.mu is not None:
            logger.info("       mu =", self.mu)
        if self.md is not None:
            logger.info("       md =", self.md)
        if self.ms is not None:
            logger.info("       ms =", self.ms)
        if self.m is not None:
            logger.info("        m =", self.m)
        if self.pre is not None:
            logger.info("      pre =", self.pre)
        if self.msml is not None:
            logger.info("    ms/ml =", self.msml)
        if self.msmu is not None:
            logger.info("    ms/mu =", self.msmu)
        if self.mdmu is not None:
            logger.info("    md/mu =", self.mdmu)
        if self.scale in ('fk', 'r1', 'r0'):
            a = self.geta()
            if self.Nt is not None:
                logger.info("        T =", round(self.getT(), 2), "[MeV]")
            logger.info("        a =", round(a, 4), "[fm]")
            if self.Ns is not None:
                logger.info("       Ls =", round(self.getLs(), 4), "1/[MeV]")
        else:
            logger.warn(f"Lattice spacing and temperature extraction not yet implemented for scale '{self.scale}'.")
        logger.info()
