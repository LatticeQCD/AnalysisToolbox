# 
# test_plot_formats.py
# 
# D. Clarke
# 
# Testing various aspects of formatting plots in the LatticeToolbox.
#

import matplotlib.pyplot as plt
from latqcdtools.base.plotting import plot_file, set_params, latexify
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

latexify()


def testPlotFormats():

    plot_file("wurf.dat", 1, 3, 4, 2, style="dots", label="wurf", marker = 'o')

    set_params(xlabel   ="$T$ [MeV]",
               ylabel   ="$\\chi^{\\rm{bare}}_{\\Re P}$",
               title    ="$N_\\tau=8$",
               legendpos=(1,0),
               legend_title="$\\ev{\\int dx}$",
               labelsintoplot=True,
               xlabelpos=(0.6,0.05),
               ylabelpos=(0.05,0.4),
               font_size=10)
    set_params(xmin=1.0,
               xmax=3.0,
               ymin=1.0
              )

    plt.show()


if __name__ == '__main__':
    testPlotFormats()