# 
# testPrintErrorBars.py                                                               
# 
# D. Clarke
# 
# Test some of the error bar stuff.
# 


from latqcdtools.base.printErrorBars import get_err_str, getValuesFromErrStr
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

def testPrintErrorBars():

    lpass = True

    x  = -5.50e-08
    xe = 2.14e-08

    if get_err_str(x,xe) != "-0.000000055(22)":
        logger.TBFail('conservative')
        lpass = False
    if get_err_str(x,xe,rounding='canonical') != "-0.000000055(21)":
        logger.TBFail('canonical')
        lpass = False

    test = "-0.000000055(22)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,-5.5e-8,ye,2.2e-8,text=test,prec=1e-8) 

    test = "-0.00000005500(22)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,-5.5e-8,ye,2.2e-10,text=test,prec=1e-8) 

    test = "3.4(1.4)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,3.4,ye,1.4,text=test,prec=1e-8) 

    test = "3.45(99)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,3.45,ye,0.99,text=test,prec=1e-8) 

    test = "1.587(32)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,1.587,ye,0.032,text=test,prec=1e-8) 

    test = "1.58007(32)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,1.58007,ye,0.00032,text=test,prec=1e-8) 

    test = "0.44(14)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,0.44,ye,0.14,text=test,prec=1e-8) 

    test = "456(14)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,456,ye,14,text=test,prec=1e-8) 

    test = "4560(140)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,4560,ye,140,text=test,prec=1e-8) 

    test = "318.45(54)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,318.45,ye,0.54,text=test,prec=1e-8) 

    test = "318.4500000(54)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,318.45,ye,0.0000054,text=test,prec=1e-8) 

    test = "-4560.(140.)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,-4560,ye,140,text=test,prec=1e-8) 

    test = "-4560(1)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,-4560,ye,1,text=test,prec=1e-8) 

    test = "0.4611(7)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,0.4611,ye,0.0007,text=test,prec=1e-8) 

    test = "4.611(7)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,4.611,ye,0.007,text=test,prec=1e-8) 

    test = "1.0(1.0)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,1.,ye,1.,text=test,prec=1e-8) 

    test = "14.0(1.0)"
    y, ye = getValuesFromErrStr(test)
    lpass *= print_results(y,14.0,ye,1.,text=test,prec=1e-8) 

    concludeTest(lpass)


if __name__ == '__main__':
    testPrintErrorBars()