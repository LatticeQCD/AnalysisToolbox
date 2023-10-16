# 
# testSpline.py                                                               
# 
# D. Clarke
# 
# Test some of the convenience wrappers for spline methods.
# 

import numpy as np
from latqcdtools.math.spline import even_knots, random_knots, getSpline
from latqcdtools.base.plotting import plt, plot_dots, plot_lines, set_params
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger
from latqcdtools.base.initialize import DEFAULTSEED


logger.set_log_level('INFO')
SHOWPLOT = False


def testSpline():

    lpass = True

    x  = np.linspace(-1, 1, 101)
    y  = 10*x**2 + np.random.randn(len(x))
    ye = np.repeat(1.,len(y))

    knots = even_knots(x, 3)

    lpass *= print_results(knots,[-0.5, 0.0, 0.5], text="even_knots")

    knots=random_knots(x, 3, SEED=DEFAULTSEED)

    lpass *= print_results(knots,[-0.56, 0.11000000000000004, 0.5000000000000001], text="random_knots")

    aicc_arr = []
    for knots in [10,30,60]:

        spline, aicc = getSpline(x,y,num_knots=knots,order=3,edata=ye,getAICc=True)

        aicc_arr.append(aicc)

        if SHOWPLOT:
            set_params(xlabel='x',ylabel='y')
            plot_dots(x, y, ye)
            plot_lines(x,spline(x),marker=None)
            plt.show()

        logger.info('Try spline with knots =',knots,'AICc =',aicc)

    if aicc_arr[0]>aicc_arr[-1]:
        logger.TBFail("AICc should penalize overfitting.")
        lpass=False
    else:
        logger.TBPass("AICc penalizes overfitting.")

    concludeTest(lpass)


if __name__ == '__main__':
    testSpline()
