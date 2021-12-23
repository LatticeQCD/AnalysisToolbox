import numpy as np
from latqcdtools.statistics import std_mean

# Ntw=Nt warm = Nt (usually)
# Ntc=Nt cold = Ntprime (usually)

def Grec(data, Ntw):
    Ntc = len(data)
    corr = [np.mean(i) for i in data]
    return Grec_direct(corr, Ntw, Ntc)


def Grec_direct(corr, Ntw, Ntc):
    res = []
    for tau in range(Ntw):
        sum = 0
        for taup in range(tau, Ntc - Ntw + tau + 1, Ntw):
            sum += corr[taup]
        res.append(sum)
    return np.array(res)



def GdivGrec_direct(g, grec):
    return g/grec


def GdivGrec(data):
    Nt_warm = len(data[1])
    grec = Grec(data[0], Nt_warm)
    g = std_mean(data[1], axis = 1)
    return GdivGrec_direct(g, grec)

def GDiffGrecDiff_direct(g, grec):
    return (g[:-1] - g[1:]) / (grec[:-1] - grec[1:])

def GDiffGrecDiff(data):
    Nt_warm = len(data[1])
    grec = Grec(data[0], Nt_warm)
    g = std_mean(data[1], axis = 1)
    return GDiffGrecDiff_direct(g, grec)

def GSubGrecSub_direct(g, grec):
    return (g - g[int(len(g)/2)]) / (grec - grec[int(len(grec)/2)])

def GSubGrecSub(data):
    Nt_warm = len(data[1])
    grec = Grec(data[0], Nt_warm)
    g = std_mean(data[1], axis = 1)
    return GSubGrecSub_direct(g, grec)


