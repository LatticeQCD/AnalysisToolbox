latqcdtools.physics.denseObs
=============

```Python
mean_square(vec):
'''Unbiased calculation of < vec**2 >. '''
```
```Python
op_to_obs(opTable, lp, obs=None, filename='denseObservables.d'):
'''
Take the operators from loadDens and combine them into physically meaningful observables. Some terminology:
    l--light
    s--strange
    B--baryon number
    Q--electric charge
    I--isospin
    S--strangeness

Parameters
----------
opTable : dict
    A table indexed by confID. Its values are numpy arrays of operators that have been measured.
lp : latticeParams
    Parameters for the ensemle the configuration belongs to.
obs : observablesOfInterest, optional
    A list of the observables you want to compute.
filename : str, optional
    Name for output table.
'''
```
```Python
class observablesOfInterest(iterable=None):
'''
A class to specify the dense observables you want to look at. It contains some consistency checks, like making
sure that you don't add an observable that is not yet computable. It also has some attributes that help streamline
using methods like np.genfromtxt. 
'''
```
