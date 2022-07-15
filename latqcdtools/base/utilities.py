#
# utilities.py
#
# D. Clarke
#
# Some utilities that you might use in any program.
#
from subprocess import run, PIPE
import time
import latqcdtools.base.logger as logger


def getArgs(parser):
    """ Get arguments from the ArgumentParser. Complain if you don't get exactly the correct arguments. """
    args, invalid_args = parser.parse_known_args()
    if len(invalid_args)>0:
        logger.TBError("Received unrecognized arguments",invalid_args)
    return args


def printArg(message,param):
    """ Some arguments are None by default, and you only want to print the if they are set. """
    if param is not None:
        print(message,param)


def shell(*args):
  """ Carry out the passed arguments args in the shell. Can be passed as a single
      string or as a list. Captures and returns output of shell command. E.g.
        shell('ls -lah')
  """
  args = [str(s) for s in args]
  process = run(' '.join(args),shell=True,check=True,stdout=PIPE,universal_newlines=True)
  return process.stdout


#
# A case where he fails:
# ['thermalTable_mu0.0357', 'thermalTable_mu0.0952', 'thermalTable_mu0.1309', 'thermalTable_mu0.0833',
#  'thermalTable_mu0.0595', 'thermalTable_mu0.0119', 'thermalTable_mu0.0', 'thermalTable_mu0.0714',
#  'thermalTable_mu0.0476', 'thermalTable_mu0.1071', 'thermalTable_mu0.119', 'thermalTable_mu0.0238']
#
def naturalSort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


class timer:

    """ A class to facilitate doing rudimentary timings in the Toolbox. """

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
