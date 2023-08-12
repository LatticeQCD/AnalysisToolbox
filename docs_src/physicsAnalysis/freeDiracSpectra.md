# Spectrum of the Dirac operator for free Wilson fermions and Mobius Domain wall fermions 
This is based on the discussions given in, 	arXiv:hep-lat/0511052. Currently, we only have Wilson fermions and Mobius Domain wall fermions. One can also come up with his own fermion discretization scheme and compare the eigenvalues.

## diracFreespectra.py

Module Contents:
1. GammaMatrix Class:
Represents the gamma matrices used in quantum field theory in Euclidean.

Methods:

g(i=1):
Input: Integer i (values between 1 to 4). Represents the specific gamma matrix desired.
Output: The corresponding 4x4 gamma matrix.
g5():
Output: Returns the 4x4 gamma5 matrix.
2. DiracOp Class (Inherits from GammaMatrix):
Represents the Dirac Operator on a spacetime lattice.

Attributes:

Lx, Ly, Lz, Lt: Lattice extents in the four spacetime directions. These should be provided upon instantiation.
fermion: The type of fermion being used. Default is "Wilson".
Methods:

p():
Output: Computes and returns the momentum values px, py, pz, pt based on the provided lattice extents.
WilsonOp(p, mass):
Input: Momentum vector p and scalar mass.
Output: Returns the computed Wilson Dirac operator.
DWMobius4D(p, mass, M=1, b=1.5, c=0.5, Ls=12):
Input: Momentum vector p, scalar mass, and optional parameters for Mobius fermions (M, b, c, Ls).
Output: Returns the Mobius Domain Wall fermion 4D Dirac operator.
eigvalues(mass, M=1, b=1.5, c=0.5, Ls=12):
Input: Scalar mass and optional parameters for Mobius fermions (M, b, c, Ls).
Output: Computes and returns the eigenvalues of the Dirac operator for the provided lattice and fermion type as an array.
Usage:
Instantiate the DiracOp class with the desired lattice dimensions and fermion type.
Use the eigvalues method to compute the eigenvalues of the Dirac operator for the lattice.
Examples:
python
Copy code
# Create a DiracOp instance for a 4x4x4x4 lattice with Wilson fermions:
D = DiracOp(Lx=4, Ly=4, Lz=4, Lt=4, fermion="Wilson")

# Compute the eigenvalues:
eigenvalues = D.eigvalues(mass=0.1)





