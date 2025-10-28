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


def mean_square(vec):
    """ Unbiased calculation of < vec**2 >. """
    N = len(vec)
    #                         diagonal terms
    return ( np.sum(vec)**2 - np.sum(vec**2) )/( N*(N-1))


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
        logger.TBRaise("This hasn't been tested for imaginary muB.")
    if lp.Nf != '21':
        logger.TBRaise("This analysis assumes Nf=2+1.")

    # Construct the output table
    for cID in opTable:

        if len(cID) != len(cID.strip()):
            logger.TBRaise('confIDs must not have whitespace! This throws off the column indexing.')

        trMdMl=np.array(opTable[cID][0]) # tr M^-1 d M 
        trMdMs=np.array(opTable[cID][1])

        if len(trMdMl) != len(trMdMs): 
            logger.warn("len(trMdMl) != len(trMdMs) cID = "+cID+"... skipping")
            continue

        trMdMl2=np.array(opTable[cID][2])  # tr ( M^-1 d M )^2
        trMdMs2=np.array(opTable[cID][3])

        if len(trMdMl2) != len(trMdMs2):
            logger.warn("len(trMdMl2) != len(trMdMs2) cID = "+cID+"... skipping")
            continue

        trMd2Ml=np.array(opTable[cID][4]) # tr ( M^-1 dd M )^2
        trMd2Ms=np.array(opTable[cID][5])

        if len(trMd2Ml) != len(trMd2Ms): 
            logger.warn("len(trMd2Ml) != len(trMd2Ms), cID = "+cID+"... skipping")
            continue

        if len(trMdMl)==0 or len(trMdMl2)==0 or len(trMd2Ml)==0:
            logger.warn("Found zero random vectors for an observable, cID = "+cID+"... skipping")
            continue

        # I follow the QCD Thermodynamics section of my researchNotes: https://github.com/clarkedavida/researchNotes
        # In the dense code, each trace comes with a 1/vol4. So whenever we have stuff like obs**2, since each factor 
        # obs has a trace, we need to multiply by vol4 to get a correct normalization. Number densities should be 
        # pure imaginary configuration by configuration at mu=0, so we take imaginary parts to reduce the noise. 
        # When this quantity is squared, it introduces a (-) sign.

        chi2l   = - vol4*( mean_square(trMdMl.imag) )/16 - np.mean(trMdMl2.real)/4 + np.mean(trMd2Ml.real)/4 + 0j
        chi2s   = - vol4*( mean_square(trMdMs.imag) )/16 - np.mean(trMdMs2.real)/4 + np.mean(trMd2Ms.real)/4 + 0j
        chi11ll = - vol4*( mean_square(trMdMl.imag) )/16 + 0*1j 
        chi11ls = - vol4*( np.mean(trMdMl.imag)*np.mean(trMdMs.imag) )/16 + 0*1j

        chi2B   = (1/9)*( 2*chi2l + chi2s + 2*chi11ll + 4*chi11ls )
        chi2Q   = (1/9)*( 5*chi2l + chi2s - 4*chi11ll - 2*chi11ls )
        chi2S   = chi2s
        chi11BQ = (1/9)*(   chi2l - chi2s +   chi11ll -   chi11ls ) 
        chi11BS = (1/3)*(         - chi2s             - 2*chi11ls ) 
        chi11QS = (1/3)*(           chi2s             -   chi11ls ) 

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