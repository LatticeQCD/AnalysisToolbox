latqcdtools.interfaces.collaborations
=============

```python
paramFromEnsLabel(ensemble, format='MILC') -> dict:
"""
Given an ensemble string, get the parameters out of it. This is meant to be paired with
the latticeParameters subclass, i.e. it gives back information to quickly construct
a latticeParameters object. 

Args:
    ensemble (str): ensemble label

Returns:
    dict: Extracted parameters 
"""
```
```python
class HotQCDParams(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf=None, scaleYear=None, muB=0):
"""
A class to handle and check the input parameters of a lattice run, especially for HotQCD.
"""
```
```python
class HotQCD_MILC_Params(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf=None, scaleYear=None, muB=0):
"""
A class to handle and check the input parameters of a lattice run using conventions common to both the
HotQCD and MILC collaborations. 
"""
```
```python
class MILCParams(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf=None, scaleYear=None, muB=0):
"""
A class to handle and check the input parameters of a lattice run, especially for MILC.
"""
```
