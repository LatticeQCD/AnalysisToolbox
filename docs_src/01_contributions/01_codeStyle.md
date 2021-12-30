# Writing Code for the ToolBox

Here are a few basic guidelines to follow when making a contribution to the ToolBox:
- Please write documentation (like this one) for any new module that you make. 
- Please write a test for each new feature you implement. We keep our tests in `AnalysisToolbox/testing`. 
- Utilize the `logger` in `latqcdtools/base/logger.py`.
- If you would like to make your python script an executable (e.g. by using `chmod +x`) you should add the line `#!/usr/bin/env python3` at the beginning of your script.
- If you are writing a method, and you want to put a comment explaining what the method generally does and how to use it, you should include it as a docstring. (You can find an explanation of docstrings "here":https://www.programiz.com/python-programming/docstrings .) An advantage to using these is that people can easily learn about your function by calling `help(functionName)`.
- If you are making a new module, please try to avoid the `from X import *` construction, as it makes it difficult to tell where function methods are truly coming from. In your main script it is of course fine. The problem with modules is just that other people may want to use them and add to them, and for that they need to understand where all the methods are defined.


# Terminal Output 

There is a dedicated logger. Its methods take arguments the same way the `print` function does. Certain methods output in colored font, for example `warn` and `TBError`. For example
```python
import latqcdtools.base.logger as logger
logger.warn("Example warning text.")
logger.TBError("Something major went wrong! This method automatically exits.")
logger.TBFail("Something went wrong. Use this for red text without exiting.")
logger.TBPass("Your test passed! Use this for green text.")
```