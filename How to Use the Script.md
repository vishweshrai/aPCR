**First read README.md**

**Primer Generator Package**

* aPCRprimer.m uses primtemp.m, primdg37.m, dGmax.m, and num2xlcol.m to generate a list of potential primers of length [18,22] nts for amplification of a target single stranded region of user-defined length and user defined circular plasmid. Choose annealing temperature based on the choice of polymerase.

* aPCR_Cycles.m estimates the final concentrations of aPCR products for the given choice of primers from the list of primers generated previously. It used the script boundtheta.m.

**Optimizer Package**

* The optimizer script apcrpareto.m runs the function apcrprimfun.m to optimize the operating conditions.
