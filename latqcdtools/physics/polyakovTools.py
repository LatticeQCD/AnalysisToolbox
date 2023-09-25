# 
# polyakovTools.py                                                               
# 
# D. Clarke
# 
# Some methods useful for calculating Polyakov loop observables. 
#

import numpy as np
from latqcdtools.base.check import checkType


# These entries come from A. Bazavov et al. "Polyakov loop in 2+1 flavor QCD from low to high temperatures".
# Phys Rev D 93 114502 (2016). Uses what they call the direct renormalization scheme. First we take the entries in
# their table, then fit that with a spline. These are the results of that.
cT={
    5.940 :    -0.38453,
    5.975 :    -0.38908,
    6.010 :    -0.39316,
    6.025 :    -0.39478,
    6.038 :    -0.39611,
    6.045 :    -0.39680,
    6.050 :    -0.39729,
    6.060 :    -0.39823,
    6.062 :    -0.39842,
    6.075 :    -0.39959,
    6.090 :    -0.40086,
    6.100 :    -0.40168,
    6.105 :    -0.40207,
    6.115 :    -0.40284,
    6.120 :    -0.40321,
    6.125 :    -0.40357,
    6.135 :    -0.40427,
    6.150 :    -0.40527,
    6.175 :    -0.40679,
    6.185 :    -0.40735,
    6.195 :    -0.40788,
    6.204 :    -0.40834,
    6.210 :    -0.40863,
    6.220 :    -0.40909,
    6.245 :    -0.41014,
    6.255 :    -0.41051,
    6.260 :    -0.41069,
    6.285 :    -0.41148,
    6.290 :    -0.41163,
    6.295 :    -0.41176,
    6.300 :    -0.41189,
    6.315 :    -0.41225,
    6.330 :    -0.41256,
    6.335 :    -0.41265,
    6.341 :    -0.41275,
    6.354 :    -0.41295,
    6.365 :    -0.41308,
    6.370 :    -0.41314,
    6.372 :    -0.41316,
    6.385 :    -0.41327,
    6.390 :    -0.41331,
    6.423 :    -0.41342,
    6.425 :    -0.41342,
    6.430 :    -0.41342,
    6.445 :    -0.41338,
    6.450 :    -0.41336,
    6.460 :    -0.41331,
    6.465 :    -0.41328,
    6.470 :    -0.41324,
    6.474 :    -0.41321,
    6.488 :    -0.41307,
    6.500 :    -0.41293,
    6.505 :    -0.41287,
    6.510 :    -0.41280,
    6.515 :    -0.41272,
    6.542 :    -0.41226,
    6.545 :    -0.41221,
    6.550 :    -0.41211,
    6.575 :    -0.41156,
    6.585 :    -0.41132,
    6.600 :    -0.41093,
    6.608 :    -0.41071,
    6.664 :    -0.40897,
    6.640 :    -0.40976,
    6.680 :    -0.40841,
    6.700 :    -0.40767,
    6.712 :    -0.40720,
    6.725 :    -0.40668,
    6.733 :    -0.40636,
    6.740 :    -0.40607,
    6.754 :    -0.40547,
    6.770 :    -0.40477,
    6.794 :    -0.40368,
    6.800 :    -0.40340,
    6.825 :    -0.40220,
    6.840 :    -0.40146,
    6.850 :    -0.40096,
    6.880 :    -0.39942,
    6.910 :    -0.39782,
    6.950 :    -0.39561,
    6.990 :    -0.39332,
    7.030 :    -0.39096,
    7.100 :    -0.38671,
    7.150 :    -0.38359,
    7.280 :    -0.37526,
}



class polyakovTools:

    """Methods for Polyakov loop observables. The methods all take as input measurements of the real and imaginary
       parts of the Polyakov loop as numpy arrays, polReIm."""

    def __init__(self, Nsigma, Ntau):
        checkType(Nsigma,int)
        checkType(Ntau,int)
        self.Ns = Nsigma
        self.Nt = Ntau


    def __repr__(self) -> str:
        return "polyakovTools"


    # Absolute value of Polyakov loop

    def absPLoop(self, polReIm):
        """ <|P|> """
        ReP=polReIm[0]
        ImP=polReIm[1]
        # Calculate |P|_ij on each configuration i belonging to the jackknife sample j, then
        # compute the mean |P|_j of that jackknife sample.
        return np.mean( np.sqrt(ReP**2+ImP**2) )


    def absPLoop2(self, polReIm):
        """ <|P|>^2 """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return np.mean( np.sqrt(ReP**2+ImP**2) )**2


    # Susceptibility and susceptibility-like observables

    def Suscept(self, polReIm):
        """ chi_|P| """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2 )


    def T3Suscept(self, polReIm):
        """ T^3 chi_|P| """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return (self.Ns/self.Nt)**3*( np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2 )


    def ReSuscept(self, polReIm):
        """ chi_ReP """
        ReP=polReIm[0]
        return self.Ns**3*( np.mean(ReP**2)-np.mean(ReP)**2 )


    def ImSuscept(self, polReIm):
        """ chi_ImP """
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ImP**2)-np.mean(ImP)**2 )


    def T3ReSuscept(self, polReIm):
        """ T^3 chi_ReP """
        ReP=polReIm[0]
        return (self.Ns/self.Nt)**3*( np.mean(ReP**2)-np.mean(ReP)**2 )


    def T3ImSuscept(self, polReIm):
        """ T^3 chi_ImP """
        ImP=polReIm[1]
        return (self.Ns/self.Nt)**3*( np.mean(ImP**2)-np.mean(ImP)**2 )


    def ReTimesIm(self, polReIm):
        """ V <ReP ImP> """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ReP*ImP) )


    def ReP2(self, polReIm):
        """ V <ReP^2> """
        ReP=polReIm[0]
        return self.Ns**3*( np.mean(ReP**2) )


    def ImP2(self, polReIm):
        """ V <ImP^2> """
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ImP**2) )


    def P2minReP2(self,polReIm):
        """ V ( <|P|>^2 - <ReP>^2 ) """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(np.sqrt(ReP**2+ImP**2))**2-np.mean(ReP)**2 )


    # Susceptibility ratios


    def RatSuscFunction(self, polReIm):
        """ R_T = chi_ImP / chi_ReP """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return (np.mean(ImP**2) - np.mean(ImP)**2) / (np.mean(ReP**2) - np.mean(ReP)**2)


    def RatSuscFunctionA(self, polReIm):
        """ R_A = chi_|P| / chi_ReP """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return ( (np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2) /
                 (np.mean(ReP**2)-np.mean(ReP)**2) )
