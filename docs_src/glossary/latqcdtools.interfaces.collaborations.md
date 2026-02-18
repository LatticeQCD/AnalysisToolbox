latqcdtools.interfaces.collaborations
=============

```Python
paramFromEnsLabel(ensemble, format='MILC'):
'''
Given an ensemble string, get the parameters out of it. 

Args:
    ensemble (str): ensemble label

Returns:
    tuple: Ns, Nt, Nf, beta string, mass1 string, mass2 string, mass3 string
'''
```
```Python
class HotQCDParams(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf='21', scaleYear=None, muB=0):
'''
A class to handle and check the input parameters of a lattice run, especially for HotQCD.
'''
```
```Python
class HotQCD_MILC_Params(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf='21', scaleYear=None, muB=0):
'''
A class to handle and check the input parameters of a lattice run using conventions common to both the
HotQCD and MILC collaborations. 
'''
```
```Python
class MILCParams(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf='21', scaleYear=None, muB=0):
'''
A class to handle and check the input parameters of a lattice run, especially for MILC.
'''
```
