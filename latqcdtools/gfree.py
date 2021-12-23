import numpy as np
from latqcdtools.num_int import int_lag_gauss, epsilon

def kernel(odT, tauT):
    return np.cosh(odT * (0.5 - tauT)) / np.sinh(odT / 2.0)

#Free Wilson correlators and spectral functions
#Source: Florian Meyer, Phd thesis

def shifted_kernel(omega, tauT, mdT):
    return np.cosh((omega + 2*mdT) * (tauT - 0.5)) / np.sinh((omega + 2*mdT) / 2.0)

def spec_pseud(omega, mdT):
    return 3 / (16.0 * np.pi**2) * (omega + 2 * mdT)**2 * np.tanh((omega + 2 * mdT)
        / 4.0) * np.sqrt(1 - ((2 * mdT)**2 / (omega + 2 * mdT)**2)) * 2

def spec_vec(omega, mdT):
    return 3 / (16.0 * np.pi**2) * (omega + 2 * mdT)**2 * np.tanh((omega + 2 * mdT)
        / 4.0) * np.sqrt(1 - ((2 * mdT)**2 / (omega + 2 * mdT)**2)) * (
            4 + 2 * (2 * mdT/(omega + 2 * mdT))**2)


def Gfree_vec(tauT, mdT):
    res = []
    for i in range(100, 180, 10):
        res.append(int_lag_gauss(lambda x: spec_vec(x, mdT)
            * shifted_kernel(x, tauT, mdT) * np.exp(x), deg = i))
    return epsilon(np.array(res))



def Gfree_pseud(tauT, mdT):
    res = []
    for i in range(100, 180, 10):
        res.append(int_lag_gauss(lambda x: spec_pseud(x, mdT)
            * shifted_kernel(x, tauT, mdT) * np.exp(x), deg = i))
    return epsilon(np.array(res))
