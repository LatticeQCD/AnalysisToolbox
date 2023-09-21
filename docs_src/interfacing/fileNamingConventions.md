# Ensemble names and file formats 

Not every collaboration names their files in the same way. On the other hand, there are a few prominent collaborations
that have their own naming schemes. It can make it easier for you to read in their files, if you can generate
these names automatically. Moreover some collaborations like to use third-party tools that have specific file
format requirements, and it is helpful to be able to read these in. Being able to interface within this highly
heterogeneous field is an important part of interoperability.

To this end, the AnalysisToolbox has some methods specifically for this purpose. At the moment we support
- MILC code configuration naming schemes
- HotQCD configuration naming schemes
- Reading `gpl` files
- Reading `yaml` files

For the first two, we have `HotQCDParams` and `MILCParams` objects, which inherit from the `latticeParams`
object described [here](../physicsAnalysis/latticeParameters.md). To quickly extract run parameters
from a MILC or HotQCD-type string, one finds inside
```Python
latqcdtools.interfaces.interfaces
```
the method `paramFrom_HotQCD_MILC`. The `gpl` and `yaml` reading methods `loadGPL` and `loadYAML`
are also inside this module.
