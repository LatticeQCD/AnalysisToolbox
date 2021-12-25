import latqcdtools.logger as logger


class O2_3d:

  " 3d O(2) critical exponents from Phys. Lett. B 492, 219 (2000). "

  d     = 3
  beta  = 0.3490
  nu    = 0.6723
  omega = 0.79
#  delta = nu*d/beta-1.
  delta = 4.7798
  alpha = 2.-beta*(1.+delta)
  gamma = beta*(delta-1.)

  def exponentSummary(self):
    print("\n Summary of 3d, O(2) critical exponents:")
    print(" alpha =",round(self.alpha,4))
    print("  beta =",round(self.beta,4))
    print(" gamma =",round(self.gamma,4))
    print(" delta =",round(self.delta,4))
    print(" omega =",round(self.omega,4))
    print("    nu =",round(self.nu,4),"\n")


class O4_3d:

  " 3d O(4) critical exponents from Nucl. Phys. B 675, 533-554 (2003). "

  d     = 3
  beta  = 0.380
  delta = 4.824
  alpha = 2.-beta*(1.+delta)
  gamma = beta*(delta-1.)
  nu    = (beta/d)*(1+delta) 
#  omega = 0.79

  def exponentSummary(self):
    print("\n Summary of 3d, O(4) critical exponents:")
    print(" alpha =",round(self.alpha,4))
    print("  beta =",round(self.beta,4))
    print(" gamma =",round(self.gamma,4))
    print(" delta =",round(self.delta,4))
#    print(" omega =",round(self.omega,4))
    print("    nu =",round(self.nu,4),"\n")


class Z2_3d:

  " 3d Z_2 critical exponents from Nucl. Phys. B 655, 277-299 (2003). "

  d     = 3
  beta  = 0.3258 
  nu    = 0.6304 
  delta = nu*d/beta-1.
  alpha = 2.-beta*(1.+delta)
  gamma = beta*(delta-1.)
#  omega = 0.79

  def exponentSummary(self):
    print("\n Summary of 3d, Z_2 critical exponents:")
    print(" alpha =",round(self.alpha,4))
    print("  beta =",round(self.beta,4))
    print(" gamma =",round(self.gamma,4))
    print(" delta =",round(self.delta,4))
#    print(" omega =",round(self.omega,4))
    print("    nu =",round(self.nu,4),"\n")


def fallFactorial(n,m):
  ''' Falling factorial n fall to m. '''
  if m==0:
    return 1
  if m>n:
    logger.TBError("m>n in falling factorial.")
  prod=1
  for i in range(m):
    prod *= n-i
  return prod


def riseFactorial(n,m):
  ''' Rising factorial n rise to m. '''
  if m==0:
    return 1
  prod=1
  for i in range(m):
    prod *= n+i
  return prod

