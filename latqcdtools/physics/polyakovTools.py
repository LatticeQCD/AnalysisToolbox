# 
# polyakovTools.py                                                               
# 
# D. Clarke
# 
# Some methods useful for calculating Polyakov loop observables. 
#

import numpy as np
from latqcdtools.base.check import checkType



class polyakovTools:

    """
    Methods for Polyakov loop observables. The methods all take as input measurements of the real and imaginary
    parts of the Polyakov loop as numpy arrays, polReIm.
    """

    def __init__(self, Nsigma, Ntau):
        checkType('int',Nsigma=Nsigma)
        checkType('int',Ntau=Ntau)
        self.Ns = Nsigma
        self.Nt = Ntau


    def __repr__(self) -> str:
        return "polyakovTools"


    # Absolute value of Polyakov loop

    def absPLoop(self, polReIm):
        """ 
        <|P|> 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        # Calculate |P|_ij on each configuration i belonging to the jackknife sample j, then
        # compute the mean |P|_j of that jackknife sample.
        return np.mean( np.sqrt(ReP**2+ImP**2) )


    def absPLoop2(self, polReIm):
        """ 
        <|P|>^2 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return np.mean( np.sqrt(ReP**2+ImP**2) )**2


    # Susceptibility and susceptibility-like observables

    def Suscept(self, polReIm):
        """ 
        chi_|P| 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2 )


    def T3Suscept(self, polReIm):
        """ 
        T^3 chi_|P| 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return (self.Ns/self.Nt)**3*( np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2 )


    def ReSuscept(self, polReIm):
        """ 
        chi_ReP 
        """
        ReP=polReIm[0]
        return self.Ns**3*( np.mean(ReP**2)-np.mean(ReP)**2 )


    def ImSuscept(self, polReIm):
        """ 
        chi_ImP 
        """
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ImP**2)-np.mean(ImP)**2 )


    def T3ReSuscept(self, polReIm):
        """ 
        T^3 chi_ReP 
        """
        ReP=polReIm[0]
        return (self.Ns/self.Nt)**3*( np.mean(ReP**2)-np.mean(ReP)**2 )


    def T3ImSuscept(self, polReIm):
        """ 
        T^3 chi_ImP 
        """
        ImP=polReIm[1]
        return (self.Ns/self.Nt)**3*( np.mean(ImP**2)-np.mean(ImP)**2 )


    def ReTimesIm(self, polReIm):
        """ 
        V <ReP ImP> 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ReP*ImP) )


    def ReP2(self, polReIm):
        """ 
        V <ReP^2> 
        """
        ReP=polReIm[0]
        return self.Ns**3*( np.mean(ReP**2) )


    def ImP2(self, polReIm):
        """ 
        V <ImP^2> 
        """
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(ImP**2) )


    def P2minReP2(self,polReIm):
        """ 
        V ( <|P|>^2 - <ReP>^2 ) 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return self.Ns**3*( np.mean(np.sqrt(ReP**2+ImP**2))**2-np.mean(ReP)**2 )


    # Susceptibility ratios


    def RatSuscFunction(self, polReIm):
        """ 
        R_T = chi_ImP / chi_ReP 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return (np.mean(ImP**2) - np.mean(ImP)**2) / (np.mean(ReP**2) - np.mean(ReP)**2)


    def RatSuscFunctionA(self, polReIm):
        """ 
        R_A = chi_|P| / chi_ReP 
        """
        ReP=polReIm[0]
        ImP=polReIm[1]
        return ( (np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2) /
                 (np.mean(ReP**2)-np.mean(ReP)**2) )
