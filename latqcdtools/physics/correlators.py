# 
# correlators.py                                                               
# 
# D. Clarke
# 
# Some tools useful when analyzing correlators 
# 
import numpy as np
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger

def foldCorrelator(corr) -> np.ndarray:
    checkType(np.ndarray,corr=corr)
    if corr.ndim != 1:
        logger.TBRaise('corr should be 1-d')
    N           = len(corr)
    half_N      = N // 2
    first_half  = corr[1:half_N]
    second_half = corr[N-1:half_N:-1]
    averaged    = (first_half + second_half) / 2.0
    result      = np.concatenate(([corr[0]], averaged, [corr[half_N]]))
    return result
