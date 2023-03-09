# Loading configurations

There are a handful of gauge configuration formats on the market, including
- [MILC](https://github.com/milc-qcd/milc_qcd)
- NERSC 
- ILDG, which you can learn about [here](https://hpc.desy.de/ildg/) and [here](https://www.sciencedirect.com/science/article/pii/S0010465511000476)

It is sometimes convenient to have the ability to read these configurations in Python. For instance lots of
machine learning code is in Python. It may also be easier to use Python to write short scripts to convert between
one format and another. To this end, we have implemented `latqcdtools/interfaces/confReader.py`.

At the moment it only reads NERSC format, but one can easily extend this, by following the example of
the `NERSCReader` class, which inherits from a more general `confReader` class. If you want to implement your own class,
please have it inherit from `confReader` as well.

A `NERSCReader` object is easily instantiated with just information about $N_s$ and $N_\tau$ as
```Python
reader = NERSCReader(Ns=8, Nt=4)
```
One then loads a gauge field using
```Python
gauge = reader.readConf('nersc.l8t4b3360')
```
You can access the link at site $(0,0,1,1)$ pointing in the 0-direction through
```Python
gauge.getLink(0,0,1,1,0)
```
