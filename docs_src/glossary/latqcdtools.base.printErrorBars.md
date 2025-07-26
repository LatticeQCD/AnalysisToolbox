latqcdtools.base.printErrorBars
=============

```Python
getValuesFromErrStr(errStr):
'''
Convert A string of the form XX.XX(YY) into a float mean and error bar. Scientific notation not yet supported.
Single digit errors not yet supported.

Args:
    errStr (str): string of the form XX.XXXX(YY).

Returns:
    float, float: mean, error. 
'''
```
```Python
get_err_str(param, param_err, numb_err_dig=2, rounding='conservative') -> str:
'''
Get the string of a number + error, e.g. 1.234567+-0.324456 --> 12.34(33) (numb_err_dig = 2). 

Args:
    param (float): _description_
    param_err (float): _description_
    numb_err_dig (int, optional): How many significant digits for your error? Defaults to 2.
    rounding (str, optional): The strategy for rounding the last significant digit of the error.
        The 'canonical' strategy rounds how one learns in grade school, i.e. 44 gets rounded to 40.
        The 'conservative' strategy always rounds up, i.e. 44 gets rounded to 50. Defaults to
        'conservative'. 

Returns:
    str: Error string. 
'''
```
```Python
get_err_str_auto(param, param_err, numb_err_dig=1, mulicon='x') -> str:
'''
Automatically express the number with an exponent. The multiplication icon can be changed. 
'''
```
```Python
get_err_str_auto_tex(param, param_err, numb_err_dig=1) -> str:
'''
Automatically express the number with an exponent in LaTeX. 
'''
```
```Python
get_err_str_exp(param, param_err, exp, numb_err_dig=1, multicon='x') -> str:
'''
Express the number with a exponent. The multiplication icon can be changed. 
'''
```
```Python
get_err_str_exp_tex(param, param_err, exp, numb_err_dig=1) -> str:
'''
Express the number with an exponent in LaTeX. 
'''
```
```Python
get_exp(param):
'''
Get the exponent of a number to base 10. 
'''
```
