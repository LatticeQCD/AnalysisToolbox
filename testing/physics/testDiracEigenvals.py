#
# testDiracEigenvals.py
#
# J. Goswami
#
# Python implementation of the spectrum of the Lattice QCD dirac operators.
#

from latqcdtools.physics.diracFreespectra import DiracOp
from latqcdtools.base.plotting import latexify, set_params, plot_dots, plt

latexify()

N = 16
Lx, Ly, Lz, Lt = N, N, N, 8
mass = 0.0


def testEigenevals():

    DW = DiracOp(Lx, Ly, Lz, Lt, fermion='Wilson')
    eig_vals = DW.eigvalues(mass)
    set_params(xlabel='${\\rm Re}\\,\\lambda$', ylabel='${\\rm Im}\\,\\lambda$')
    plot_dots(eig_vals.real, eig_vals.imag)
    plt.show()

    MDW = DiracOp(Lx, Ly, Lz, Lt, fermion="DwMobius")
    eig_vals = MDW.eigvalues(mass)
    set_params(xlabel='${\\rm Re}\\,\\lambda$', ylabel='${\\rm Im}\\,\\lambda$')
    plot_dots(eig_vals.real, eig_vals.imag)
    plt.show()


if __name__ == '__main__':
    testEigenevals()