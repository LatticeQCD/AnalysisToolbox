Contributions
=============

Here are a few basic guidelines to follow when making a contribution to the AnalysisToolbox:
- Please **write documentation** (like this one) for any new module that you make. 
- Please **write a test for each new feature you implement**. We keep our tests in `AnalysisToolbox/testing`.
- Make sure you **run all the tests** when you're done to check that nothing was broken. You can find more information
about the tests [here](testing.md). 
- For all screen output and for terminating your program because of errors: Please try to use the `logger` in `latqcdtools/base/logger.py`. 
This has the advantage that all scripts are terminated in the same way. Also you can have your program make little log files in 
addition to the terminal output, in case you want that! Finally in the case of error, warning, and debug messages, the logger will
report to you where the messages are coming from.
- For all plotting: Please try to use `latqcdtools/base/plotting.py`.
- For reading or writing any data table: Please try to use `readTable` and `writeTable` in `latqcdtools/base/readWrite.py`.
- If you would like to make your python script an executable (e.g. by using `chmod +x`) you should add the line `#!/usr/bin/env python3` 
at the beginning of your script.
- Please **be mindful of backward compatibility**. David is sometimes the worst offender in this regard. In particular, try not to
delete functions, even if it seems like they are not being used. If you find a better way of organizing things, you may need
to implement wrappers to maintain backward compatibility. There are many people using the AnalysisToolbox, and we don't want
to disrupt their workflows too much.
- If you are writing a method, and you want to put a comment explaining what the method generally does and how to use it, you should include 
it as a *docstring*. (You can find an explanation of docstrings [here](https://www.programiz.com/python-programming/docstrings).) 
An advantage to using these is that people can easily learn about your function by calling `help(functionName)`. Certain IDEs, such as 
VSCode, include tools to generate your docstrings largely automatically. Also in VSCode, when I hover over a function, it reports
to me its docstring, which I find helpful.
- **We do not allow the `from X import *` construction** anywhere within the AnalysisToolbox, as it makes it difficult to tell where 
function methods are truly coming from. The problem with modules is just that other people may want to use them and add to them, 
and for that they need to understand where all the methods are defined. This syntax also confuses IDEs.
- Please wrap all your `main_*` applications and `test*` tests like so:
  ```Python
  def myMainOrTest():
     ...
  
  if __name__ == '__main__':
      myMainOrTest()
  ```
  This has a few advantages: (1) you can put functions outside this definition, then import your main program and borrow those 
  functions without running anything; (2) the logger will give more meaningful output when called in your main function; and (3) 
  weirdly this seems to help some of our parallelizers run on some Macs.
- Please **include a `__repr__` function for all your classes**. For instance
  ```Python
  class CAT():
        def __init__(self):
            """ A class for cats. """
            ...
        def __repr__(self) -> str:
            return "CAT"
  ```
  The `__repr__` function is a special method that provides a string representation of an instance of the class. This is used by
  for instance the [logger](logging.md), which when printing e.g. `DEBUG` messages, will report what method (or class if applicable)
  called it. 
- The `-> str` syntax in the above example indicates that we expect a string as a return value. I don't think this is used by the Python
  interpreter, but it helps keep the code self-documenting. You don't have to use such *return type annotations*, because I think in
  general you may write methods that might return multiple kinds of types. In the case of `__repr__` there is only one possibility,
  so I try to include it.
- If you have a method that you think shouldn't be called outside of a module, for instance some kind of internal helper method
  or auxiliary method, try to prefix that method with an underscore, like `_method()`.
- Please for the love of God **don't use lambda functions**. I've been coding in Python since about 2012, and I 
  still can't read these things. Why can't you just implement it as a normal function? Using a lambda function saves one line, 
  but drastically reduces readability. 

Here are some more detailed articles:

```{toctree}
---
maxdepth: 1
---
testing.md
git.md
documentation.md
logging.md
```
