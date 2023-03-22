Contributions
=============

Here are a few basic guidelines to follow when making a contribution to the AnalysisToolBox:
- Please write documentation (like this one) for any new module that you make. 
- Please write a test for each new feature you implement. We keep our tests in `AnalysisToolbox/testing`. 
- For all screen output and for terminating your program because of errors: Please try to use the `logger` in `latqcdtools/base/logger.py`. This has the advantage that all scripts are terminated in the same way. Also you can have your program make little log files in addition to the terminal output, in case you want that! 
- For all plotting: Please try to use `latqcdtools/base/plotting.py`.
- For reading or writing any data table: Please try to use `readTable` and `writeTable` in `latqcdtools/base/readWrite.py`.
- If you would like to make your python script an executable (e.g. by using `chmod +x`) you should add the line `#!/usr/bin/env python3` at the beginning of your script.
- If you are writing a method, and you want to put a comment explaining what the method generally does and how to use it, you should include it as a docstring. (You can find an explanation of docstrings [here](https://www.programiz.com/python-programming/docstrings).) An advantage to using these is that people can easily learn about your function by calling `help(functionName)`.
- We do not allow the `from X import *` construction anywhere within the AnalysisToolbox, as it makes it difficult to tell where function methods are truly coming from. The problem with modules is just that other people may want to use them and add to them, and for that they need to understand where all the methods are defined.
- Please wrap all your `main_*` applications and `test*` tests like so:
  ```Python
  def myMainOrTest():
     ...
  
  if __name__ == '__main__':
      myMainOrTest()
  ```
  This has a few advantages: (1) you can put functions outside this definition, then import your main program and borrow those functions without running anything; (2) the logger will give more meaningful output when called in your main function; and (3) weirdly this seems to help things run on some Macs.

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
