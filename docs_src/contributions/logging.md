# Logging


There is a dedicated logger in
```Python
import latqcdtools.base.logger as logger
```
Its methods take arguments the same way the `print` function does. Certain methods output in colored font, for 
example `warn` and `TBError`. For instance
```Python
logger.warn("Example warning text.")
logger.TBError("Something major went wrong! This method automatically exits.")
logger.TBFail("Something went wrong. Use this for shiny red text without exiting.")
logger.TBPass("Your test passed! Use this for shiny green text.")
```

Most of our methods use this logger. Some advantages of using it are
- Automatic timestamps
- Optional output to log file, which works with multiple processors
- Recording the git commit

If you would like to save your output to a log file, we recommend using the module
```Python
latqcdtools.base.initialize
```
You can start a log file for your run using
```Python
initialize('myRun.log')
```
at the beginning of your code.
In this case, all the output will get print to screen as well as to `myRun.log`.
This includes the git hash by default, which is useful for debugging and reproducability.
You can also call `finalize()` at the end of your code, but for now no real cleanup is required,
so this just prints an "I'm finished!" message.

