#
# timer.py
#
# D. Clarke
#
import time


class timer:

    """A class to facilitate doing rudimentary timings in the Toolbox."""

    def __init__(self):
        print("\n  Timer initialized.\n")
        self._tstart = time.time()
        self._tend   = self._tstart

    def printTiming(self):
        self._tstart = self._tend
        self._tend   = time.time()
        timing = self._tend - self._tstart
        print("  Time to finish: %12.8f [s]" % timing)
