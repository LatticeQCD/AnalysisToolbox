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
various information criteria, Bayesian model averaging, 
estimation of integrated autocorrelation time,
error ellipses, Kolmogorov-Smirnov tests, 
and curve fitting with and without Bayesian priors.
We stress that our math and statistics methods are generally useful, independent of physics contexts.
- **General physics:** Unit conversions, critical exponents and temperatures for various universality
classes, physical constants, framework for spin models. 
- **QCD physics:** Hadron resonance gas model, ideal fermion gas, HotQCD equation of state, 
and the QCD beta function. 
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


To use the AnalysisToolbox, make sure you have Python 3.9+. You should then be able to
conveniently install it using
```bash
pip install latqcdtools
```
If you would like alternative ways to set up the AnalysisToolbox, please see
[the documentation](https://latticeqcd.github.io/AnalysisToolbox).

Besides this basic install, there is a `latexify()` command you can use when plotting to 
make your plot font match typical LaTeX documents. In order for this command to work, you need
to have LaTeX installed on your system. The easiest is to install `texlive-full`, but
if that is not possible, it may be enough to install `texlive-mathscience` in addition
to the basic stuff.


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
[K. Ebira](https://github.com/kaiebira),
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
- DAC acknowledges helpful discussions with 
[C. DeTar](https://github.com/detar),
[X.-Y. Jin](https://github.com/jxy),
[S. Lahert](https://github.com/Shaun252), 
[G. P. LePage](https://github.com/gplepage), 
and [C. Rohde](https://github.com/rohdog2003).

