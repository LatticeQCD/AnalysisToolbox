latqcdtools.interfaces.HotQCD
=============

```Python
loadDens(densFile, confID, lp, inTable=None) -> dict:
'''
Allows reading of output from C. Schmidt's Dense code. The Dense code produces measurements of various operators
relevant for calculating conserved charge fluctuations. We store as a dictionary indexed by confID, which lets us
conveniently combine data when there are multiple dense files per configuration. Here we update the table
inTable for the new read-in, which yields outTable. Only supports Nf=2+1 configurations for the time being.

Parameters
----------
densFile : str
    Name of the Dense code file to be read.
confID : str
    A unique identifier for the configuration on which dense calculation was performed. Every person seems to have a
    different scheme for labelling configurations, so this needs to be a string to be as flexible as possible.
lp : latticeParams
    Parameters for the ensemle the configuration belongs to.
inTable : dict
    A table indexed by confID. Its values are a list of operators that have been measured.

Returns
-------
outTable : dict
    A table indexed by confID. Its values are a list of operators that have been measured.
'''
```
```Python
makeConfTag(conf, stream) -> str:
'''
This takes a configuration number conf and stream label stream to make a tag labelling a configuration.
Implementing this as a function makes sure everyone using the Toolbox as the same convention and, more importantly,
ensures that the tags have no whitespace in them, which otherwise can throw off the column counting of methods in
the denseObs module. 
'''
```
```Python
class HotQCDParams(Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=None, Nf='21', scaleYear=None, muB=0):
'''
A class to handle and check the input parameters of a lattice run, especially for HotQCD.
'''
```
