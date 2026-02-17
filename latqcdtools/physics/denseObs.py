# 
# denseObs.py                                                               
# 
# D. Clarke
# 
# Methods to turn output from C. Schmidt's dense code into observables.
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import writeTable
from latqcdtools.interfaces.HotQCD import getObs

def mean_square(vec):
    """ 
    Unbiased calculation of < vec**2 >. 
    """
    N = len(vec)
    #                         diagonal terms
    return ( np.sum(vec)**2 - np.sum(vec**2) )/( N*(N-1))

def square_mean_unequal(vec1,vec2):
    """
    <vec1><vec2> when vec1 and vec2 are of unequal length. This is needed because some
    vectors receive extra statistics.
    """
    N1 = len(vec1)
    N2 = len(vec2)
    N = np.min([N1,N2])
    return np.mean(vec1[:N])*np.mean(vec2[:N])

def op_to_obs(opTable,lp,writeFiles=True,outFolder='denseObservables') -> dict:
    """
    Take the operators from loadDens and combine them into physically meaningful observables. Some terminology:
        l--light
        s--strange
        B--baryon number
        Q--electric charge
        I--isospin
        S--strangeness

    Args:
        opTable (dict): Operators loaded from loadDens. 
        lp (HotQCD_MILC_Params): latticeParams object with info about the ensemble. 
        writeFiles (bool, optional): Write final observables in denseObservables directory. Defaults to True.

    Returns:
        dict: Final observables 
    """

    OBS = {
    "confID" : [],
#    "nl"     : [],  # net light-quark density
#    "ns"     : [],  # net strange-quark density
#    "nB"     : [],  # net baryon-number density
#    "nQ"     : [],  # net electric-charge density
#    "nS"     : [],  # net strangeness-number density
#    "nI"     : [],  # net isospin density
#    "nl^2"   : [], 
#    "nB^2"   : [],  
#    "nQ^2"   : [], 
#    "nS^2"   : [],
    "X2l"    : [],  # chi_2^l
    "X2s"    : [],  # chi_2^s
    "X11ll"  : [],  # chi_11^ll
    "X11ls"  : [],  # chi_11^ls
    "X2B"    : [],  # chi_2^B
    "X2Q"    : [],  # chi_2^Q
    "X2S"    : [],  # chi_2^S
    "X11BQ"  : [],
    "X11BS"  : [],
    "X11QS"  : [],
    }

    vol4 = lp.vol4
    if lp.muB != 0:
        logger.TBRaise("This analysis assumes muB=0.") 
    if lp.Nf != '21':
        logger.TBRaise("This analysis assumes Nf=2+1.")

    # Construct the output table
    for cID in opTable:

        try:
            trMdMl     , trMdMs      = getObs(opTable,cID,'trMdM')
            trMdMl2    , trMdMs2     = getObs(opTable,cID,'trMdM2')
            trMd2Ml    , trMd2Ms     = getObs(opTable,cID,'trMd2M')
            trMd3Ml    , trMd3Ms     = getObs(opTable,cID,'trMd3M')
            trMd4Ml    , trMd4Ms     = getObs(opTable,cID,'trMd4M')
            trMd2M2l   , trMd2M2s    = getObs(opTable,cID,'trMd2M2')
            trMdMMd2Ml , trMdMMd2Ms  = getObs(opTable,cID,'trMdMMd2M')
            trMdMMd3Ml , trMdMMd3Ms  = getObs(opTable,cID,'trMdMMd3M')
            trMdM3l    , trMdM3s     = getObs(opTable,cID,'trMdM3')
            trMdM4l    , trMdM4s     = getObs(opTable,cID,'trMdM4')
            trMdM2Md2Ml, trMdM2Md2Ms = getObs(opTable,cID,'trMdM2Md2M')
        except logger.ToolboxException:
            continue

        # Express in terms of A operators. mu enters with an i in the partition function, so when I take derivatives
        # w.r.t. mu, it pulls down an i factor. When I then evaluate at mu=0, that i is all that remains. Hence odd
        # derivatives are pure imaginary while even derivatives are pure real. Moreover, when I square pure imaginary
        # operators, it becomes real again, with a minus sign. 

        trAl     = trMdMl.imag     # second order
        trAs     = trMdMs.imag
        trA2l    = trMd2Ml.real
        trA2s    = trMd2Ms.real
        trAe2l   = trMdMl2.real
        trAe2s   = trMdMs2.real

        trAe3l   = trMdM3l.imag    # third order
        trAe3s   = trMdM3s.imag
        trAA2l   = trMdMMd2Ml.imag
        trAA2s   = trMdMMd2Ms.imag
        trA3l    = trMd3Ml.imag # is this zero?
        trA3s    = trMd3Ms.imag

        trAe4l   = trMdM4l.real    # fourth order 
        trAe4s   = trMdM4s.real
        trAe2A2l = trMdM2Md2Ml.real 
        trAe2A2s = trMdM2Md2Ms.real
        trA2e2l  = trMd2M2l.real
        trA2e2s  = trMd2M2s.real
        trAA3l   = trMdMMd3Ml.real 
        trAA3s   = trMdMMd3Ms.real
        trA4l    = trMd4Ml.real 
        trA4s    = trMd4Ms.real 

        trdAl    = trA2l - trAe2l
        trdAs    = trA2s - trAe2s
        trddAl   = trA3l - trAA2l - 2*trAe3l
        trddAs   = trA3s - trAA2s - 2*trAe3s
        trddA2l  = trA4l - 2*trAA3l - trA2e2l + 2*trAe2A2l
        trddA2s  = trA4s - 2*trAA3s - trA2e2s + 2*trAe2A2s

        # I follow the QCD Thermodynamics section of my researchNotes: https://github.com/clarkedavida/researchNotes
        # In the dense code, each trace comes with a 1/vol4. So whenever we have stuff like obs**2, since each factor 
        # obs has a trace, we need to multiply by vol4 to get a correct normalization. 
        chi2l   = - vol4*( mean_square(trAl) )/16 - np.mean(trAe2l)/4 + np.mean(trA2l)/4 + 0j
        chi2s   = - vol4*( mean_square(trAs) )/16 - np.mean(trAe2s)/4 + np.mean(trA2s)/4 + 0j
        chi11ll = - vol4*( mean_square(trAl) )/16 + 0j 
        chi11ls = - vol4*( np.mean(trAl)*np.mean(trAs) )/16 + 0j

        chi2B   = (1/9)*( 2*chi2l + chi2s + 2*chi11ll + 4*chi11ls )
        chi2Q   = (1/9)*( 5*chi2l + chi2s - 4*chi11ll - 2*chi11ls )
        chi2S   = chi2s
        chi11BQ = (1/9)*(   chi2l - chi2s +   chi11ll -   chi11ls ) 
        chi11BS = (1/3)*(         - chi2s             - 2*chi11ls ) 
        chi11QS = (1/3)*(           chi2s             -   chi11ls ) 

        # I guess something here is wrong for fourth order 
        chi4l = (1/2)*np.mean(trAe4l + 7*trAe2A2l - 2*trA2e2l - 3*trAA3l + trA4l) \
                -(1/8)*(np.mean(trddAl)*np.mean(trAl)+np.mean(trdAl)) \
                +(1/8)*(square_mean_unequal(trddAl,trAl) + mean_square(trdAl))

# Idea on why imaginary part is wrong: in your eq (12.22), the <> brackets are the gauge average. The tr is already
# a random vector average, and that random vector average is exactly what you're computing in here. You don't do the
# gauge average until the very end, way outside of this routine, so when you needed to compute the <>**2 term, you
# accidentally computed <()**2>, which is not the same obviously. This still works at real chemical potential
# because that term is zero, so you never needed to include it.
#
# The way around this would be to identify all the pieces that enter an observable, like mean_square(trAl), then
# give those configuration by configuration. Then one can correctly reconstruct the jackknife averages also for
# pure imaginary chemical potential.

# TODO: Find some data to test these against. They should be correct, but
#       it's better to be careful.
#        nl2  = - mean_square(trMdMl.imag)*vol4/16 + 0j
#        ns2  = - mean_square(trMdMs.imag)*vol4/16 + 0j
#        nlns = np.mean(trMdMl.imag)*np.mean(trMdMs.imag)*vol4/16 +0j
#        nl   =  np.mean( trMdMl )/4
#        ns   =  np.mean( trMdMs )/4
#        nB   =  ( 2*nl  + ns           )/3
#        nB2  =  ( 4*nl2 + ns2 + 4*nlns )/9
#        nl2  =      nl2
#        dnl  = -( np.mean(trMdMl2) - np.mean(trMd2Ml) )/4
#        dns  = -( np.mean(trMdMs2) - np.mean(trMd2Ms) )/4
#        dnS  =           dns
#        dnQ  = ( 5*dnl + dns )/9
#        dnI  =     dnl + dns
#        dnB  = ( 2*dnl + dns )/9
#        nS   = -        ns
#        nQ   =  (  nl - ns)/3
#        nI   = complex(0)
#        nS2  =            ns2
#        nQ2  =  (   nl2 + ns2 - 2*nlns )/9

        OBS["confID"].append(    cID         )
#        OBS["nl"    ].append(     nl*lp.Nt   )
#        OBS["ns"    ].append(     ns*lp.Nt   )
#        OBS["nB"    ].append(     nB*lp.Nt   )
#        OBS["nQ"    ].append(     nQ*lp.Nt)
#        OBS["nS"    ].append(     nS*lp.Nt)
#        OBS["nI"    ].append(     nI*lp.Nt) 
#        OBS["nl^2"  ].append(    nl2*lp.Nt**2)
#        OBS["nB^2"  ].append(    nB2*lp.Nt**2) 
#        OBS["nQ^2"  ].append(    nQ2*lp.Nt**2)
#        OBS["nS^2"  ].append(    nS2*lp.Nt**2)
        OBS["X2l"   ].append(  chi2l*lp.Nt**2)
        OBS["X2s"   ].append(  chi2s*lp.Nt**2)
        OBS["X11ll" ].append(chi11ll*lp.Nt**2)
        OBS["X11ls" ].append(chi11ls*lp.Nt**2)
        OBS["X2B"   ].append(  chi2B*lp.Nt**2)
        OBS["X2Q"   ].append(  chi2Q*lp.Nt**2)
        OBS["X2S"   ].append(  chi2S*lp.Nt**2)
        OBS["X11BQ" ].append(chi11BQ*lp.Nt**2)
        OBS["X11BS" ].append(chi11BS*lp.Nt**2)
        OBS["X11QS" ].append(chi11QS*lp.Nt**2)

    if writeFiles:
        logger.info(f"Write observables in {outFolder}/{lp.getcparams()}...")
        for observable in OBS:
            if len(OBS[observable])>0 and observable!="confID":
                logger.info("  ",observable)
                writeTable(f'{outFolder}/{lp.getcparams()}/{observable}.txt',OBS["confID"],OBS[observable],header=["confID",observable])

    return OBS