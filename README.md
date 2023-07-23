# LatticeToolbox

The LatticeToolbox is a collection of Python tools that can be used for data analysis, with an aim in 
particular toward lattice QCD. The LatticeToolbox is used especially by the HotQCD collaboration, so
it also includes some methods complying with their conventions. It was originally created by
H. Sandmeyer during his PhD work, and has since been maintained by some of the HotQCD scientists.

## Getting Started

In order to use these scripts, please ensure that you have the following:
1. Python 3.6 or better
2. colorama
3. cycler
4. matplotlib
5. numba
6. numpy
7. pathos
8. scipy
9. sympy
10. LaTeX (Probably best if you install TeXLive-Full)


For your convenience, packages (2-9) can be installed via
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

## Documentation

Please check out [the documentation](https://latticeqcd.github.io/LatticeToolbox) to learn how to use the 
LatticeToolbox.


## Getting help and bug report
Open an [issue](https://github.com/LatticeQCD/LatticeToolbox/issues), if...
- you have troubles running the code.
- you have questions on how to implement your own routine.
- you have found a bug.
- you have a feature request.

If none of the above cases apply, you may also send an email to clarke(dot)davida(at)gmail(dot)com.


## Contributors

[D. Clarke](https://github.com/clarkedavida), 
[L. Altenkort](https://github.com/luhuhis), 
[J. Goswami](https://github.com/jishnuxx),
[L. Mazur](https://github.com/lukas-mazur),
[H. Sandmeyer](https://github.com/hsandmeyer),
[M. Sarkar](https://github.com/mugdhasarkar),
[C. Schmidt](https://github.com/schmidt74), 
[H.-T. Shu](https://github.com/haitaoshu), 
T. Ueding

## Acknowledgment

- We acknowledge support by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) through the CRC-TR 211 'Strong-interaction matter under extreme conditions'– project number 315477589 – TRR 211.
- This work was partly performed in the framework of the PUNCH4NFDI consortium supported by DFG fund "NFDI 39/1", Germany.
- This work is also supported by the U.S. Department of Energy, Office of Science, though the Scientific Discovery through Advance Computing (SciDAC) award
Computing the Properties of Matter with Leadership Computing Resources.

