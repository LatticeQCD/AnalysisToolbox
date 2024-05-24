# AnalysisToolbox

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://latticeqcd.github.io/AnalysisToolbox)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/LatticeQCD/AnalysisToolbox/commits/main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8259640.svg)](https://doi.org/10.5281/zenodo.8241876)



The AnalysisToolbox set of Python tools for statistically analyzing correlated data. 
This includes aspects of lattice QCD applications related to QCD phenomenology. 


We advertise briefly here some features of the AnalysisToolbox:
- **General mathematics:** Numerical differentiation, convenience wrappers for
SciPy numerical integration and solving IVPs.
- **General statistics:** Jackknife, bootstrap, Gaussian bootstrap, error propagation,
various information criteria, estimation of integrated autocorrelation time,
error ellipses, Kolmogorov-Smirnov tests, 
and curve fitting with and without Bayesian priors.
We stress that our math and statistics methods are generally useful, independent of physics contexts.
- **General physics:** Unit conversions, critical exponents for various universality
classes, physical constants, framework for spin models. 
- **QCD physics:** Hadron resonance gas model, HotQCD equation of state, and the QCD beta function. 
These methods are useful for QCD phenomenology, independent of lattice contexts.
- **Lattice QCD:** Continuum-limit extrapolation, Polyakov loop observables, SU(3) gauge
fields, reading in gauge fields, and the static quark-antiquark potential. These methods
rather target lattice QCD. 


In any of the above cases, after installing the AnalysisToolbox, you can easily incorporate
its features in your own Python scripts like any other library.
Some simple examples are in the [tutorial](https://latticeqcd.github.io/AnalysisToolbox/tutorial.html).
A realistic use-case that weaves the AnalysisToolbox into a lattice
QCD workflow can be found in [this data publication](https://pub.uni-bielefeld.de/record/2979080).
More information can be found in [the documentation](https://latticeqcd.github.io/AnalysisToolbox).

## Setting up the AnalysisToolbox

In order to use these scripts, please ensure that you have the following:
1. Python 3.9+
2. colorama
3. cycler
4. matplotlib
5. numba
6. numpy
7. pathos
8. pyyaml
9. scipy
10. sympy
11. pytest
12. LaTeX (Probably best if you install TeXLive-Full)

You can set up the AnalysisToolbox using either a [Python virtual environment](https://realpython.com/python-virtual-environments-a-primer/)
or by modifying the Python on your local machine. You can set up the former or the latter by adjusting the user preferences
in `configureToolbox.py`.

### Configuring for a virtual environment

The recommended way to proceed is to create a Python virtual environment 
in which you install all the required packages. This is what the Python people seem to prefer, which you can
read about in their [PEP 668](https://peps.python.org/pep-0668/). It is also more convenient for use on
supercomputers, since you don't need superuser privileges to install any missing packages.
On the other hand, virtual environments can be slower.

In `configureToolbox.py` set
```Python
STRATEGY = "VENV" 
```
then run
```Bash
./configureToolbox.py
```
This will create a `venv` folder containing all the information about your virtual environment. Every time you open a new
terminal, if you want to use the AnalysisToolbox, you will need to 
```Bash
cd scripts
source startVENV.bash
```
You can terminate your virtual environment any time using
```Bash
deactivate
```

### Configuring using your OS Python directly 

If you're old-fashioned like David is, you can also just directly `pip3 install` on your machine,
modifying your OS Python.

In `configureToolbox.py` set
```Python
STRATEGY = "BASIC" 
```
then run
```Bash
./configureToolbox.bash
```

### Installing the required packages

Either strategy will make sure your `PYTHONPATH` environment variable points
to the correct place. You will need to close your terminal and open a new one.

Once you carried out one of the above two strategies,
packages (2-10) can be installed via
```shell
pip3 install -r requirements.txt
```
There are some further packages required if you would like to make contributions to the AnalysisToolbox; in particular
there are many packages needed to compile the documentation. If you are interested in writing documentation, you should also 
instead
```shell
pip3 install -r docRequirements.txt
```

Once this has all been settled, try running the tests using
```shell
pytest
```


## Getting started and documentation

To acquaint yourself with the AnalysisToolbox, you can start by
having a look at the [tutorial](https://latticeqcd.github.io/AnalysisToolbox/tutorial.html),
which walks through some scripts in the `examples` directory.
You can also look at some of the scripts in the `applications` and `testing` directories.

To learn about the code in more detail, especially learning how to contribute, please have
a look [the documentation](https://latticeqcd.github.io/AnalysisToolbox).


## Getting help and bug reports

Open an [issue](https://github.com/LatticeQCD/AnalysisToolbox/issues), if...
- you have troubles running the code.
- you have questions on how to implement your own routine.
- you have found a bug.
- you have a feature request.

If none of the above cases apply, you may also send an email to clarke(dot)davida(at)gmail(dot)com.


## Contributors

[D. A. Clarke](https://github.com/clarkedavida), 
[L. Altenkort](https://github.com/luhuhis), 
[H. Dick](https://github.com/redweasel),
[J. Goswami](https://github.com/jishnuxx),
[O. Kaczmarek](https://github.com/olaf-kaczmarek),
[L. Mazur](https://github.com/lukas-mazur),
[H. Sandmeyer](https://github.com/hsandmeyer),
[M. Sarkar](https://github.com/mugdhasarkar),
[C. Schmidt](https://github.com/schmidt74), 
[H.-T. Shu](https://github.com/haitaoshu), 
[T. Ueding](https://github.com/SiggiUphues)

## Crediting AnalysisToolbox

If you used this code in your research, your teaching, or found it generally useful, please help
us out by citing
```
@inproceedings{Altenkort:2023xxi,
    author = "Altenkort, Luis and Clarke, David Anthony and Goswami, Jishnu and Sandmeyer, Hauke",
    title = "{Streamlined data analysis in Python}",
    eprint = "2308.06652",
    archivePrefix = "arXiv",
    primaryClass = "hep-lat",
    month = "8",
    year = "2023"
}
```

## Acknowledgments

- We acknowledge support by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) through the CRC-TR 211 'Strong-interaction matter under extreme conditions'– project number 315477589 – TRR 211.
- This work was partly performed in the framework of the PUNCH4NFDI consortium supported by DFG fund "NFDI 39/1", Proj.No. 460248186 (PUNCH4NFDI).
- DAC acknowledges helpful discussions with [C. DeTar](https://github.com/detar), S. Lahert, 
and [G. P. LePage](https://github.com/gplepage).

