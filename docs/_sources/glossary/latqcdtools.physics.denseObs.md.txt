latqcdtools.physics.denseObs
=============

```Python
mean_square(vec):
'''Unbiased calculation of < vec**2 >. '''
```
```Python
op_to_obs(opTable, lp, writeFiles=True, outFolder='denseObservables') -> dict:
'''
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
'''
```
