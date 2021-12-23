import numpy as np
import latqcdtools.solve as sv


def tmp(meff, nt, Nt, cn, cnp):
    return corr_cosh(meff,nt,Nt)-cn/cnp

def corr_cosh(meff, nt, Nt):
    return np.cosh(meff * (nt - Nt / 2)) / np.cosh(meff * (nt + 1 - Nt / 2))


def eff_cosh_mass(cn, cnp, nt, Nt, x1 = 0, x2 = 5, tol = 1e-14):
    try:
        return sv.solve(corr_cosh, cn / cnp, x1, x2, tol, nt, Nt)
    except ValueError:
        return np.nan


def calc_eff_mass_direct(mean, xdata, Nt):
    res = []
    for i, nt in enumerate(xdata[:-1]):
        try:
            if mean[i] / mean[i + 1] < 0 or (mean[i] / mean[i + 1] > 1
                    and i >= Nt / 2):
                res.append(np.nan)
            else:
                effm = eff_cosh_mass(mean[i], mean[i + 1], nt, Nt)
                if effm == 0:
                    res.append(np.nan)
                else:
                    res.append(effm)
        except IndexError:
            res.append(np.nan)
        
    return np.array(res)


def calc_eff_mass(data, xdata, Nt):
    mean = [np.mean(data[i]) for i in range(len(data))]
    return calc_eff_mass_direct(mean, xdata, Nt)
