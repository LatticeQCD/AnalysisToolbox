# AnalysisToolbox

The AnalysisToolbox is a collection of Python tools that can be used for data analysis, with an aim in particular toward lattice QCD.

## Prerequisites

In order to use these scripts, please ensure that you have the following:
1. Python 3.6
2. numpy
3. scipy
4. mpmath
5. matplotlib
6. LaTeX
7. colorama (for colors in logging messages)
8. sklearn
9. tqdm (for making progress bars)

Additionally you need to have python3 somewhere in your $PATH. To do this, you have to define the environment variable PYTHONPATH containing the path to the root folder of this project; e.g. add the line<pre><code class="shell">export PYTHONPATH="${PYTHONPATH}:/path/to/your/AnalysisToolbox/"
</code></pre>to your .bashrc.
