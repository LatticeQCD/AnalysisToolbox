# 
# polyakovTools.py                                                               
# 
# D. Clarke, 2 Mar 2021 
# 
# Some methods useful for calculating Polyakov loop observables. 
#
import numpy as np

class polyakovTools:

  """Methods for Polyakov loop observables. The methods expect numpy arrays."""

  def __init__(self, Nsigma, Ntau):
    self.Ns = Nsigma
    self.Nt = Ntau 


  # Absolute value of Polyakov loop

  ''' <|P|> '''
  def absPLoop(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    # Calculate |P|_ij on each configuration i belonging to the jackknife sample j, then
    # compute the mean |P|_j of that jackknife sample.
    return np.mean( np.sqrt(ReP**2+ImP**2) )


  ''' <|P|>^2 '''
  def absPLoop2(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return np.mean( np.sqrt(ReP**2+ImP**2) )**2


  # Susceptibility and susceptibility-like observables  

  ''' chi_|P| '''
  def Suscept(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return self.Ns**3*( np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2 )


  ''' T^3 chi_|P| '''
  def T3Suscept(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return (self.Ns/self.Nt)**3*( np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2 )


  ''' chi_ReP '''
  def ReSuscept(self, polReIm):
    ReP=polReIm[0]
    return self.Ns**3*( np.mean(ReP**2)-np.mean(ReP)**2 )


  ''' chi_ImP '''
  def ImSuscept(self, polReIm):
    ImP=polReIm[1]
    return self.Ns**3*( np.mean(ImP**2)-np.mean(ImP)**2 )


  ''' T^3 chi_ReP '''
  def T3ReSuscept(self, polReIm):
    ReP=polReIm[0]
    return (self.Ns/self.Nt)**3*( np.mean(ReP**2)-np.mean(ReP)**2 )


  ''' T^3 chi_ImP '''
  def T3ImSuscept(self, polReIm):
    ImP=polReIm[1]
    return (self.Ns/self.Nt)**3*( np.mean(ImP**2)-np.mean(ImP)**2 )


  ''' V <ReP ImP> '''
  def ReTimesIm(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return self.Ns**3*( np.mean(ReP*ImP) )


  ''' V <ReP^2> '''
  def Re2(self, polReIm):
    ReP=polReIm[0]
    return self.Ns**3*( np.mean(ReP**2) )


  ''' V <ImP^2> '''
  def Im2(self, polReIm):
    ImP=polReIm[0]
    return self.Ns**3*( np.mean(ImP**2) )


  ''' V ( <|P|>^2 - <ReP>^2 ) '''
  def P2minReP2(self,polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return self.Ns**3*( np.mean(np.sqrt(ReP**2+ImP**2))**2-np.mean(ReP)**2 )


  # Susceptibility ratios  


  ''' R_T = chi_ImP / chi_ReP '''
  def RatSuscFunction(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return (np.mean(ImP**2) - np.mean(ImP)**2) / (np.mean(ReP**2) - np.mean(ReP)**2)


  ''' R_A = chi_|P| / chi_ReP '''
  def RatSuscFunctionA(self, polReIm):
    ReP=polReIm[0]
    ImP=polReIm[1]
    return ( (np.mean(ReP**2+ImP**2)-np.mean(np.sqrt(ReP**2+ImP**2))**2) /
             (np.mean(ReP**2)-np.mean(ReP)**2) )

