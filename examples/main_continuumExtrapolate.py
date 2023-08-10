import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.initialize import initialize, finalize
from latqcdtools.math.num_deriv import diff_deriv
from latqcdtools.math.spline import getSpline
from latqcdtools.statistics.statistics import gaudif
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.physics.continuumExtrap import continuumExtrapolate
from latqcdtools.physics.constants import r0_phys
from latqcdtools.physics.lattice_params import latticeParams

initialize('cont.log')

Nts         = [6,8,10,12,14,16,18,20]
Tds, Tderrs = [], []

for Nt in Nts:

    T  = []
    Ns = Nt*3

    # Read in Polyakov loop measurements, 
    beta, PM, PE = readTable('ploop/Nt'+str(Nt)+'.txt',usecols=(0,1,2))

    # Create array of temperatures in physical units
    for b in beta:
        lp = latticeParams(Ns, Nt, b, scaleType='r0') 
        T.append( lp.getT() )
    t = np.linspace(T[0],T[-1],1001)

    # Extract Td from inflection point of <|P|> vs T using natural spline
    def getTd(pm):
        spl  = getSpline(T, pm, natural=True)
        dPdT = diff_deriv(t, spl)
        maxIndex = np.argmax(dPdT)
        return t[maxIndex]

    # Error in Td estimate comes from 1000 Gaussian bootstrap samples
    Td, Tde = bootstr_from_gauss(getTd, PM, PE, 1000)
    Tds.append(Td)
    Tderrs.append(Tde)

# Perform O(a^4) continuum-limit extrapolation
result, result_err, chidof = continuumExtrapolate( Nts, Tds, Tderrs, order=2, xtype="Nt",
                                                   show_results=True, plot_results=True )

# Do a Z-test against literature result,  
r0 = r0_phys(year=2014, units="MeVinv")
Tdr0, Tdr0e = r0 * result[0], r0 * result_err[0]
Tdr0_lit, Tdr0_lite = 0.7457, 0.0045
logger.info('q(ours vs. lit) =',gaudif(Tdr0,Tdr0e,Tdr0_lit,Tdr0_lite))

finalize()
