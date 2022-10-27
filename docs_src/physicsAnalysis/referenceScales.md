# Reference scales

In this module you find lots of functions for scale setting, along with references for each scale determination. 
Supported scales include: $f_K$, $t_0$, $r_0$, $r_1$, and $r_1m_s$. 

You can get scales in lattice units
through methods like `a_times_fk_2014(beta)`, which takes the bare coupling as input and returns $a\,f_K(\beta)$
using a parameterization from a 2014 paper.
Please note that the parameterizations for $a$ have closed intervals $[\beta_{\rm min},~\beta_{\rm max}]$ in
which they are valid.

For each scale we have methods giving values in physical units along with the determination year. For example you
can get the $r_1$ scale in [fm] using `r1_MILC_2010('fm')` along with its corresponding combined error bar
`r1err_MILC_2010('fm')`.

Have a look at this module to see what other methods are available.
