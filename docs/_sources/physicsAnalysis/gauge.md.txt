# Gaugefields 

The `gaugeField` object, found in `latqcdtools/physics/gauge.py`, is a collection of `SU(3)` objects, 
described [here](../math/SU3.md). This is useful
whenever one needs to store a configuration in memory. Such objects can be instantiated as
```Python
gauge = gaugeField(self.Ns,self.Nt)
```
Its links can be read and manipulated through the accessors `getLink` and `setLink`. The `gaugeField` object is the
place to set up basic observables. For instance one can already compute the plaquette using
```Python
gauge.getPlaquette()
```
