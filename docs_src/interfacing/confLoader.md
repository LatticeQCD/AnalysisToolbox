# Loading configurations

There are a handful of gauge configuration formats on the market, including
- [MILC](https://github.com/milc-qcd/milc_qcd)
- NERSC 
- ILDG, which you can learn about [here](https://hpc.desy.de/ildg/) and [here](https://www.sciencedirect.com/science/article/pii/S0010465511000476)

It is sometimes convenient to have the ability to read these configurations in Python. For instance lots of
machine learning code is in Python. It may also be easier to use Python to write short scripts to convert between
one format and another. To this end, we have implemented `latqcdtools/interfaces/confReader.py`.

At the moment it reads NERSC and ILDG format, as well as the formats written by the dyniso code, but one can easily extend this, by following the example of e.g.
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

If you do not know the format or size of a configuration, `detectConf` finds them from the file:
```Python
fmt, Ns, Nt = detectConf('nersc.l8t4b3360')   # ('NERSC', 8, 4)
```

## Writing configurations

`latqcdtools/interfaces/confWriter.py` writes a `gaugeField` in NERSC format:
```Python
writer = NERSCWriter(gauge)
writer.writeConf('out.nersc')
```
By default the writer stores all 3 rows of each link in big-endian double precision; the `format`, `endianness`, and
`precision` arguments change this. The checksum in the header is computed from the link values as in
[SIMULATeQCD](https://github.com/LatticeQCD/SIMULATeQCD), so it does not depend on the endianness. The application
`main_confConvert.py` combines `detectConf`, the readers, and the writer to convert between formats.

## dyniso

The `DynisoReader` reads the three formats of dyniso's `save_lattice`: the compressed format, in which dyniso also
writes its checkpoints, the MILC format, and the ascii format. You do not need to say which one; the reader detects
it from the file:
```Python
reader = DynisoReader(Ns=8, Nt=4)
gauge = reader.readConf('chk0')
```
The compressed format stores each link $U$ through the eight numbers $a_k = \mathrm{Im}\,\mathrm{tr}(\lambda_k \log U)$,
with $\lambda_k$ the Gell-Mann matrices, and the reader rebuilds $U = \exp(\tfrac{i}{2}\sum_k a_k\lambda_k)$. dyniso
writes all three formats in single precision, and the ascii format with 6 significant digits, so the links are only
that precise. Like dyniso, the reader reunitarizes the links after reading them.
