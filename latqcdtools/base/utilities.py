#
# utilities.py
#
# D. Clarke
#
# Some utilities that you might use in any program.
#
import time
import latqcdtools.base.logger as logger


def getArgs(parser):
    """Get arguments from the ArgumentParser. Complain if you don't get exactly the correct arguments."""
    args, invalid_args = parser.parse_known_args()
    if len(invalid_args)>0:
        logger.TBError("Received unrecognized arguments",invalid_args)
    return args


class timer:

    """A class to facilitate doing rudimentary timings in the Toolbox."""

    def __init__(self):
        print("\n  Timer initialized.\n")
        self._tstart = time.time()
        self._tend   = self._tstart


    def printTiming(self, message=None):
        self._tstart = self._tend
        self._tend   = time.time()
        timing = self._tend - self._tstart
        if message is None:
            print("\n  Time to finish: %12.8f [s]\n" % timing)
        else:
            print("\n  "+message+"\n")


    def resetTimer(self):
        self._tstart = time.time()
        self._tend   = self._tstart
