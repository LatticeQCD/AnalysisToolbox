# LatticeToolbox

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://latticeqcd.github.io/LatticeToolbox)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/LatticeQCD/LatticeToolbox/commits/main)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7994982.svg)](https://doi.org/10.5281/zenodo.7994982) -->

The LatticeToolbox is a collection of Python tools that can be used for data analysis, with an aim in 
particular toward lattice QCD.


## Setting up the LatticeToolbox

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
11. LaTeX (Probably best if you install TeXLive-Full)

The recommended way to proceed is to create a [Python virtual environment](https://realpython.com/python-virtual-environments-a-primer/),
in which you install all the required packages. This is what the Python people seem to prefer, which you can
read about in their [PEP 668](https://peps.python.org/pep-0668/).
If you're old-fashioned like David is, you can also just directly `pip3 install` on your machine,
potentially modifying your OS Python.

Either way, for your convenience, packages (2-10) can be installed via
```shell
pip3 install -r requirements.txt
```
There are some further packages required if you would like to make contributions to the LatticeToolbox; in particular
there are many packages needed to compile the documentation. If you are interested in helping develop, you should
instead
```shell
pip3 install -r developerRequirements.txt
```
Additionally you need to make sure your `PYTHONPATH` environment variable points
to the correct place. This can be accomplished by
```shell
python3 installToolbox.py
```
You then need to close your terminal and open a new one.
Once this has all been settled,
try running the tests. You can do this by going to the `testing` folder
and calling
```shell
bash runTests.bash
```

## Getting started and documentation

To acquaint yourself with the LatticeToolbox, you can start by
having a look at the [tutorial](https://latticeqcd.github.io/LatticeToolbox/tutorial.html),
which walks through some scripts in the `examples` directory.
You can also look at some of the scripts in the `applications` and `testing` directories.

To learn about the code in more detail, especially learning how to contribute, please have
a look [the documentation](https://latticeqcd.github.io/LatticeToolbox).


## Getting help and bug reports

Open an [issue](https://github.com/LatticeQCD/LatticeToolbox/issues), if...
- you have troubles running the code.
- you have questions on how to implement your own routine.
- you have found a bug.
- you have a feature request.

If none of the above cases apply, you may also send an email to clarke(dot)davida(at)gmail(dot)com.


## Contributors

[D. A. Clarke](https://github.com/clarkedavida), 
[L. Altenkort](https://github.com/luhuhis), 
[J. Goswami](https://github.com/jishnuxx),
[O. Kaczmarek](https://github.com/olaf-kaczmarek),
[L. Mazur](https://github.com/lukas-mazur),
[H. Sandmeyer](https://github.com/hsandmeyer),
[M. Sarkar](https://github.com/mugdhasarkar),
[C. Schmidt](https://github.com/schmidt74), 
[H.-T. Shu](https://github.com/haitaoshu), 
[T. Ueding](https://github.com/SiggiUphues)

## Crediting LatticeToolbox

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
- This work was partly performed in the framework of the PUNCH4NFDI consortium supported by DFG fund "NFDI 39/1", Germany.

