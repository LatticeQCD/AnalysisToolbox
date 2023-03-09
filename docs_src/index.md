AnalysisToolbox
===========

The AnalysisToolbox is a collection of Python tools written for statistical data analysis, with an aim specifically 
toward lattice field theory. It contains a variety of useful methods for these purposes, for example general jackknife 
and bootstrap routines that can calculate error bars of arbitrary functions of data. It also includes modules tailored 
for lattice field theory, such as methods for scale setting and universal quantities from commonly examined 
3-d universality classes.

The source code is inside the `latqcdtools` folder. It is organized into several subfolders:
- `applications`: Ready-to-use applications.
   Discussed in more detail in the [applications](applications/applications.md) section.
- `base`: Fundamental classes and methods that are common to many other methods. 
   Discussed in more detail in the [base](base/base.md) section.
- `experimental`: Modules that were used in an older version of the AnalysisToolbox that are not yet part of this new
   version. These may be incorporated eventually.
- `interfaces`: Allows interfacing with other entities present in lattice QCD. This includes things like
   collaboration-specific naming conventions for configurations, as well as output into niche data file formats.
   The interfacing methods may require installation of packages that are not needed for the rest of the toolbox.
- `math`: Modules and classes for doing general math without any physics context.
- `physics`: Physics-focused modules and classes. This includes the HRG classes, as well as classes that make analyzing
   lattice quantities like the Polyakov loop easier. Discussed in more detail in
   Discussed in more detail in the [dataAnalysis](dataAnalysis/dataAnalysis.md) section.
- `statistics`: Modules and classes for doing general statistics without any physics context.
   Discussed in more detail in the [physicsAnalysis](physicsAnalysis/physicsAnalysis.md) section.
- `scripts`: Bash scripts that help write comments for AnalysisToolbox code, or to help repair it.
- `testing`: Python scripts that test the AnalysisToolbox methods.

We would love it if you are interested in helping develop the AnalysisToolbox! Please have a look to the
[contributions](contributions/contributions.md) section to learn how to do this in a nice way.

```{toctree}
---
maxdepth: 1
---
contributions/contributions.md
base/base.md
math/math.md
dataAnalysis/dataAnalysis.md
physicsAnalysis/physicsAnalysis.md
interfacing/interfacing.md
applications/applications.md
```
