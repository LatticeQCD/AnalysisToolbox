#
# testDiracEigenvals.py
#
# J. Goswami
#
# Python implementation of the spectrum of the Lattice QCD dirac operators.
#

from latqcdtools.physics.diracFreespectra import DiracOp
from latqcdtools.base.plotting import set_params, plot_dots, plt
import latqcdtools.base.logger as logger


N = 16
Lx, Ly, Lz, Lt = N, N, N, 8
mass = 0.0

SHOWPLOTS = False
if not SHOWPLOTS:
    logger.warn('Suppressing plt.show() for pytest purposes.')


def testEigenevals():

    DW = DiracOp(Lx, Ly, Lz, Lt, fermion='Wilson')
    eig_vals = DW.eigvalues(mass).ravel()
    set_params(xlabel='${\\rm Re}\\,\\lambda$', ylabel='${\\rm Im}\\,\\lambda$')
    plot_dots(xdata=eig_vals.real, ydata=eig_vals.imag)
    if SHOWPLOTS:
        plt.show()

    MDW = DiracOp(Lx, Ly, Lz, Lt, fermion="DwMobius")
    eig_vals = MDW.eigvalues(mass).ravel()
    set_params(xlabel='${\\rm Re}\\,\\lambda$', ylabel='${\\rm Im}\\,\\lambda$')
    plot_dots(xdata=eig_vals.real, ydata=eig_vals.imag)
    if SHOWPLOTS:
        plt.show()


if __name__ == '__main__':
    testEigenevals()
