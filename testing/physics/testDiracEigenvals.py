#
# testDiracEigenvals.py
#
# J. Goswami
#
# Python implementation of the spectrum of the Lattice QCD dirac operators.
#

import pylab as plt
from latqcdtools.physics.diraFreespectra import DiracOp
from multiprocessing import Pool

N = 16
Lx, Ly, Lz, Lt = N, N, N, 16
mass = 0.25

def Teseignevals():
    DW = DiracOp(Lx, Ly, Lz, Lt, fermion='Wilson')
    eig_vals = DW.eigvalues(mass)
    plt.plot(eig_vals.real, eig_vals.imag, 'r.', mfc='None')
    plt.show()

    MDW = DiracOp(Lx, Ly, Lz, Lt, fermion="DwMobius")
    eig_vals = MDW.eigvalues(mass)
    plt.plot(eig_vals.real, eig_vals.imag, 'r.', mfc='None')
    plt.show()

if __name__ == '__main__':
    Teseignevals()