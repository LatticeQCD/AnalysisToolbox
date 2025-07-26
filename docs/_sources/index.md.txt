AnalysisToolbox
===========

The AnalysisToolbox is a collection of Python tools written for statistical data analysis, with an aim specifically 
toward lattice field theory. It contains a variety of useful methods for these purposes, for example general jackknife 
and bootstrap routines that can calculate error bars of arbitrary functions of data. It also includes modules tailored 
for lattice field theory, such as methods for scale setting and universal quantities from commonly examined 
3-d universality classes. 

Ready-to-use applications are in the `applications` directory, discussed in more detail 
in the [applications](applications/applications.md) section.

The source code is inside the `latqcdtools` folder. It is organized into several subfolders:
- `base`: Fundamental classes and methods that are common to many other methods. 
   Discussed in more detail in the [base](base/base.md) section.
- `interfaces`: Allows interfacing with other entities present in lattice QCD. This includes things like
   collaboration-specific naming conventions for configurations, as well as output into niche data file formats.
   The interfacing methods may require installation of packages that are not needed for the rest of the toolbox.
- `math`: Modules and classes for doing general math without any physics context.
- `physics`: Physics-focused modules and classes. This includes the HRG classes, as well as classes that make analyzing
   lattice quantities like the Polyakov loop easier. Discussed in more detail in
   Discussed in more detail in the [dataAnalysis](dataAnalysis/dataAnalysis.md) section.
- `statistics`: Modules and classes for doing general statistics without any physics context.
   Discussed in more detail in the [physicsAnalysis](physicsAnalysis/physicsAnalysis.md) section.

Also at the highest level is a `scripts` directory, containing Bash scripts that help write comments for 
AnalysisToolbox code, or to help repair it. The `examples` directory contains some pedagogical examples
how to use the AnalysisToolbox. Finally there is a `tests` folder, which has 
unit tests for the AnalysisToolbox methods. 

We would love it if you are interested in helping develop the AnalysisToolbox! Please have a look to the
[contributions](contributions/contributions.md) section to learn how to do this in a nice way.


# Setting up the AnalysisToolbox

To use the AnalysisToolbox, make sure you have Python 3.9+. You should then be able to
conveniently install it using
```bash
pip install latqcdtools
```
Besides this, there is a `latexify()` command you can use when plotting to make your
plot font match typical LaTeX documents. In order for this command to work, you need
to have LaTeX installed on your system. We recommend installing TeXLive-Full.

## What if pip install fails?

If the `pip install` doesn't work for you, for example because you don't have access to Python 3.9,
you should be able to get something to work by cloning the AnalysisToolbox git.
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
python3 configureToolbox.py
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


```{toctree}
---
maxdepth: 1
---
tutorial.md
contributions/contributions.md
base/base.md
math/math.md
dataAnalysis/dataAnalysis.md
physicsAnalysis/physicsAnalysis.md
interfacing/interfacing.md
applications/applications.md
glossary/glossary.md
```
