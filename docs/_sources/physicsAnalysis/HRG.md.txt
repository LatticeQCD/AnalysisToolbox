# Hadron resonance gas

The hadron resonance gas is a non-interacting model, valid below the crossover temperature, that uses hadrons and
their resonances as the only degrees of freedom. In this model, the pressure for a single species $i$ comes out to be
$$
  P_i=\frac{m^2g_iT^2}{2\pi^2}\sum_{k=1}^\infty\frac{\eta_i^{k+1}z_i^k}{k^2}
                         K_2\left(\frac{m_ik}{T}\right),
$$
where $g_i$ is the species' degeneracy, $m_i$ its mass, $\eta_i$ is +1 if it is a boson and -1 if it is a fermion,
$K_2$ is the modified Bessel function of the second kind,
$$
  z_i=\exp{\mu_B B_i + \mu_Q Q_i + \mu_S S_i + \mu_C C_i},
$$
and $B_i$, $Q_i$, $S_i$, and $C_i$ are the species' baryon number, electric charge, strangeness, and charm, 
respectively. The total pressure is then
$$
  \frac{P}=\sum_{i=1}P_i,
$$
where the summation runs over states.
From this pressure, all other thermodynamic observables can be derived analytically, then evaluated
numerically

## HRG.py

In this module, we have the `HRG` class and its associated methods. This class can be instantiated using some information
from some list of hadrons and resonances, for example those lists given in the `HRGtables` subfolder.
In this class are implemented a handful of dimensionless observables and some of their temperature derivatives.

In evaluating the above sum over $k$, we use the fact that the modified Bessel functions are exponentially suppressed
in mass and only keep at most 20 terms for numerical efficiency. For states with somewhat high masses (e.g. baryons and
heavier states) we may only keep the first term. This is the so-called Boltzmann approximation.
