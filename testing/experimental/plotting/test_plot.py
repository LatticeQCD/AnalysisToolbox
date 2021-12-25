#!/usr/bin/env python3

import numpy as np
from latqcdtools.readin import *
from latqcdtools.plotting import *


latexify()
#set_default_param(alpha=0.5)

plot_file("data1.txt", 1, 2, 3, 4, label="gauss noise", xlabel="$\\frac{1}{z}$", capsize = 2)
plot_file("wurf.dat", 1, 3, 4, 2, style="lines", label="quadratic", marker = 'o')
plot_file("data2.txt", 1, 2, 3, 4, style="fill", alpha=0.3, label="linear")
plot_file("data3.txt", 1, 2, 3, 4, style="fill", alpha=0.9, label="linear 2", legendpos = "best", capsize = 2)
plot_bar([-6,-5],[1,2])
plot_func(np.sin)
plot_file_zoom('40%', '30%', 2., 4., "wurf.dat", 1, 3, 4, 2, capsize = 2, style = 'lines')
plt.savefig("tmp.pdf")
plt.close()
print("\nPlotting done. Compare tmp.pdf with ref.pdf")
