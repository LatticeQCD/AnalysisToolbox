import sys

_python34 = (sys.version_info.major == 3 and sys.version_info.minor >= 4)
_min_numpy_version = '1.11.0'
_min_scipy_version = '0.17.0'
_min_matplot_version = '1.5.1'


if not _python34:
    raise ImportError("latqcdtools requires python 3.4 or later")

try:
    import numpy
except ImportError:
    raise ImportError("latqcdtools requires numpy %s or later" % (_min_numpy_version,))

try:
    import scipy
except ImportError:
    raise ImportError("latqcdtools requires scipy %s or later" % (_min_scipy_version,))

try:
    import matplotlib
except ImportError:
    raise ImportError("latqcdtools requires matplotlib %s or later" % (_min_matplot_version,))

