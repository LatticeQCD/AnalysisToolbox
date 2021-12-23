#!/usr/bin/env python3
from latqcdtools.logger import *
from test_print import *

log_level=['DEBUG','DETAILS','PROGRESS','INFO','WARN','NONE']

for level in log_level:
    set_log_level(level)
    print('The log level is set to: ',level )
    test_print()

TBPass("Here's an","example","pass message.",[1,2,3],1.2,42)
TBError("Here's an example error message.",[1,2,3],1.2,42)
