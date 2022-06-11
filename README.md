# AnalysisToolbox

The AnalysisToolbox is a collection of Python tools that can be used for data analysis, with an aim in 
particular toward lattice QCD. The Analysistoolbox is used especially by the HotQCD collaboration, so
it also includes some methods complying with their conventions.

## Prerequisites

In order to use these scripts, please ensure that you have the following:
1. Python 3.6 or better
2. numpy
3. scipy
4. mpmath
5. matplotlib
6. colorama 
7. cycler
8. numba
9. LaTeX


For your convenience, these packages can be installed via
```shell
pip3 install -r requirements.txt
```
There are some further packages required if you would like to make contributions to the AnalysisToolBox; in particular
there are many packages needed to compile the documentation. If you are interested in helping develop, you should
instead
```shell
pip3 install -r developerRequirements.txt
```

Additionally you need to have python3 somewhere in your $PATH. To do this, you have to define the environment variable 
PYTHONPATH containing the path to the root folder of this project; e.g. add the 
```shell
export PYTHONPATH="${PYTHONPATH}:/path/to/your/AnalysisToolbox/"
```
to your `bashrc`.


## Documentation

Please check out [the documentation](https://latticeqcd.github.io/AnalysisToolbox) to learn how to use the 
AnalysisToolbox.


## Getting help and bug report
Open an [issue](https://github.com/LatticeQCD/AnalysisToolbox/issues), if...
- you have troubles running the code.
- you have questions on how to implement your own routine.
- you have found a bug.
- you have a feature request.

If none of the above cases apply, you may also send an email to dclarke(at)physik(dot)uni-bielefeld(dot)de.


## Contributors

[H. Sandmeyer](https://github.com/hsandmeyer), 
[L. Altenkort](https://github.com/luhuhis), 
[D. Clarke](https://github.com/clarkedavida), 
[J. Goswami](https://github.com/jishnuxx),
[L. Mazur](https://github.com/lukas-mazur), 
[C. Schmidt](https://github.com/schmidt74), 
[M. Sarkar](https://github.com/mugdhasarkar),
[H.-T. Shu](https://github.com/haitaoshu), 
T. Ueding

## Acknowledgment

- We acknowledge support by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) through the CRC-TR 211 'Strong-interaction matter under extreme conditions'– project number 315477589 – TRR 211.
- This work was partly performed in the framework of the PUNCH4NFDI consortium supported by DFG fund "NFDI 39/1", Germany.
- This work is also supported by the U.S. Department of Energy, Office of Science, though the Scientific Discovery through Advance Computing (SciDAC) award
Computing the Properties of Matter with Leadership Computing Resources.

