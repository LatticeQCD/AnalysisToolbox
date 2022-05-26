# 
# statisticalPhysics.py                                                               
# 
# D. Clarke
# 
# A collection of methods relevant for statistical physics calculations.
#
import latqcdtools.base.logger as logger


def printExponent(prefix, exponent):
    if exponent is not None:
        print(prefix, round(exponent, 4))


class UniversalityClass:

    """ Skeleton universality class from which all others inherit."""

    symm  = None
    d     = None
    alpha = None
    beta  = None
    gamma = None
    delta = None
    nu    = None
    eta   = None
    omega = None

    def name(self):
        return str(self.d)+"d, "+str(self.symm)

    def exponentSummary(self):
        print("\n Summary of "+self.name()+" critical exponents:\n")
        printExponent(" alpha =",self.alpha)
        printExponent("  beta =",self.beta)
        printExponent(" gamma =",self.gamma)
        printExponent(" delta =",self.delta)
        printExponent(" omega =",self.omega)
        printExponent("    nu =",self.nu)
        print()

    def hyperscalingCheck(self, tol=1e-6):
        if (2*self.beta+self.gamma-2+self.alpha) > tol:
            logger.TBError(self.name(),"fails hyperscaling check.")
        if (2*self.beta*self.delta-self.gamma-2+self.alpha) > tol:
            logger.TBError(self.name(),"fails hyperscaling check.")
        if (self.nu*self.d-2+self.alpha) > tol:
            logger.TBError(self.name(),"fails hyperscaling check.")


class O2_3d(UniversalityClass):
    """ 3d O(2) critical exponents from Phys. Lett. B 492, 219 (2000). """
    symm  = "O(2)"
    d     = 3
    beta  = 0.3490
    nu    = 0.6723
    omega = 0.79
    delta = 4.7798
    alpha = 2.-beta*(1. + delta)
    gamma = beta*(delta-1.)


class O4_3d(UniversalityClass):
    """ 3d O(4) critical exponents from Nucl. Phys. B 675, 533-554 (2003). """
    symm  = "O(4)"
    d     = 3
    beta  = 0.380
    delta = 4.824
    alpha = 2.-beta*(1.+delta)
    gamma = beta*(delta-1.)
    nu    = (beta/d)*(1+delta)


class Z2_3d(UniversalityClass):
    """ 3d Z_2 critical exponents from Nucl. Phys. B 655, 277-299 (2003). """
    symm  = "Z_2"
    d     = 3
    beta  = 0.3258
    nu    = 0.6304
    delta = nu*d/beta-1.
    alpha = 2.-beta*(1.+delta)
    gamma = beta*(delta-1.)


class Z2_2d(UniversalityClass):
    """ Exact solution for 2d Z_2 class. """
    symm  = "Z_2"
    d     = 2
    alpha = 0
    beta  = 1/8
    gamma = 7/4
    delta = 15
    nu    = 1
    eta   = 1/4