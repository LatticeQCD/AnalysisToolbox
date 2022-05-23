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
