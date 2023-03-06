# 
# testPrintErrorBars.py                                                               
# 
# D. Clarke
# 
# Test some of the error bar stuff.
# 


from latqcdtools.base.printErrorBars import get_err_str
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

def testPrintErrorBars():

    x  = -5.50e-08
    xe = 2.16e-08

    if get_err_str(x,xe) == "-0.000000055(22)":
        logger.TBPass('get_err_str')
    else:
        logger.TBFail('get_err_str')

if __name__ == '__main__':
    testPrintErrorBars()