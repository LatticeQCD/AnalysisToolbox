latqcdtools.physics.staticPotential
=============

`V_Teq0(r) -> float`

Zero temperature quark potential in [MeV], takes r in [fm]. The parameters a, b, and c come from
a Levenberg-Marquardt fit of the data in Fig 14 of Phys. Rev. D90 (2014) 094503. These numbers can
be obtained again by running analysistoolbox/hisq_potential/fit_hisq_pot.py. 

`_cpu_impdist(Ns, r2max, improvedAction=True)`

CPU implementation of improved distances calculation.
Used as fallback when CUDA is not available.

Args:
    Ns: Spatial extension of lattice
    r2max: Maximum squared distance to improve
    improvedAction: Whether to use improved action
    
Returns:
    List of improved distances

`fitV_Teq0(r, a, b, c) -> float`

Fit form of standard Cornell potential. Fit to be done in lattice units.

`fitV_Teq0_oneloop(r, a, b, c, d) -> float`

Including one-loop corrections to Coulomb. See Nucl. Phys. B 129 (1977) and Phys. Lett B 92 (1980).
Fit to be done in lattice units.

`fitV_Teq0_twoloop(r, a, b, c, d, e) -> float`

Including two-loop corrections to Coulomb. See Nucl. Phys. B 501 (1997). Fit to be done in lattice units.

`get_optimal_block_size()`

Returns an optimal block size based on the current CUDA device.
Defaults to 256 if device information cannot be obtained.

Returns:
    int: Optimal threads per block

`impdist(Ns, r2max, improvedAction=True)`

GPU-accelerated calculation of tree-level improved distances.

Follows equation (3) of 10.1103/PhysRevD.90.074038.
Falls back to CPU implementation if CUDA is unavailable.

Args:
    Ns: Spatial extension of lattice
    r2max: Maximum squared distance to improve
    improvedAction: Whether to use improved action (default: True)

Returns:
    rimp: List of improved distances

Raises:
    ValueError: If Ns <= 0 or r2max is too large

