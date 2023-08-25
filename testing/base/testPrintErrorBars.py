# 
# testPrintErrorBars.py                                                               
# 
# D. Clarke
# 
# Test some of the error bar stuff.
# 


from latqcdtools.base.printErrorBars import get_err_str, getValuesFromErrStr
from latqcdtools.math.math import print_results
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

def testPrintErrorBars():

    x  = -5.50e-08
    xe = 2.16e-08

    if get_err_str(x,xe) != "-0.000000055(22)":
        logger.TBFail('get_err_str')

    test = "-0.000000055(22)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,-5.5e-8,ye,2.2e-8,text=test,prec=1e-8) 

    test = "3.4(1.4)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,3.4,ye,1.4,text=test,prec=1e-8) 

    test = "3.45(99)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,3.45,ye,0.99,text=test,prec=1e-8) 

    test = "0.44(14)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,0.44,ye,0.14,text=test,prec=1e-8) 

    test = "456(14)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,456,ye,14,text=test,prec=1e-8) 

    test = "4560(140)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,4560,ye,140,text=test,prec=1e-8) 

    test = "318.45(54)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,318.45,ye,0.54,text=test,prec=1e-8) 

    test = "-4560.(140.)"
    y, ye = getValuesFromErrStr(test)
    print_results(y,-4560,ye,140,text=test,prec=1e-8) 


if __name__ == '__main__':
    testPrintErrorBars()