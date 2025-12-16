#
# testDiracEigenvals.py
#
# J. Goswami
#
# Python implementation of the spectrum of the Lattice QCD dirac operators.
#

import pytest
from latqcdtools.physics.diracFreespectra import DiracOp
from latqcdtools.base.plotting import set_params, plot_dots, plt


N = 16
Lx, Ly, Lz, Lt = N, N, N, 8
mass = 0.0


#@pytest.mark.skip(reason="Jishnu please fix the test, then remove this decorator when you're done.")
def testEigenevals():

    DW = DiracOp(Lx, Ly, Lz, Lt, fermion='Wilson')
    eig_vals = DW.eigvalues(mass).ravel()
    set_params(xlabel='${\\rm Re}\\,\\lambda$', ylabel='${\\rm Im}\\,\\lambda$')
    plot_dots(xdata=eig_vals.real, ydata=eig_vals.imag)
    plt.show()

    MDW = DiracOp(Lx, Ly, Lz, Lt, fermion="DwMobius")
    eig_vals = MDW.eigvalues(mass).ravel()
    set_params(xlabel='${\\rm Re}\\,\\lambda$', ylabel='${\\rm Im}\\,\\lambda$')
    plot_dots(xdata=eig_vals.real, ydata=eig_vals.imag)
    plt.show()


if __name__ == '__main__':
    testEigenevals()
