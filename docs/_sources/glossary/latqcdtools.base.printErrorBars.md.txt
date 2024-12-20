latqcdtools.base.printErrorBars
=============

`getValuesFromErrStr(errStr)`
 
    Convert A string of the form XX.XX(YY) into a float mean and error bar. Scientific notation not yet supported.
    Single digit errors not yet supported.

    Args:
        errStr (str): string of the form XX.XXXX(YY).

    Returns:
        float, float: mean, error. 
    
`get_err_str(param, param_err, numb_err_dig=2) -> str`
 
    Get the string of a number + error, e.g. 1.234567+-0.324456 --> 12.34(33) (numb_err_dig = 2). 
    
`get_err_str_auto(param, param_err, numb_err_dig=1, mulicon='x') -> str`
 
    Automatically express the number with an exponent. The multiplication icon can be changed. 
    
`get_err_str_auto_tex(param, param_err, numb_err_dig=1) -> str`
 
    Automatically express the number with an exponent in LaTeX. 
    
`get_err_str_exp(param, param_err, exp, numb_err_dig=1, multicon='x') -> str`
 
    Express the number with a exponent. The multiplication icon can be changed. 
    
`get_err_str_exp_tex(param, param_err, exp, numb_err_dig=1) -> str`
 
    Express the number with an exponent in LaTeX. 
    
`get_exp(param)`
 
    Get the exponent of a number to base 10. 
    
