# asmr-copoly
The repository contains the code needed to replicate the results of the (under review) article:

Analysis of Copolymerization with Simultaneous Reversibility and Transesterification using Stochastic Model Regression. Kuehster, L.; Dai, J.; Thompson, A.; Jhon, Y. K.; Wang, Y.; Qin, B.; Smith, W. C.; Xu, X.; Zhang, F.; Lynd, N. A. Macromolecules 2024.

Files included in this package: fit_stochastic_reversible_di_v1.cpp –
stochastic model regression for PLGA
fit_stochastic_reversible_di_mo_v1.cpp – stochastic model regression for
PLCL polymer.cpp – supporting functions polymer.h – header file for
polymer.cpp inputPLGA7525.in – example input file for PLGA
outputPLGA75-25.out – example output file for PLGA

How to use these codes: These codes enact stochastic model regression
for either PLGA or PLCL. Since lactide and glycolide are both cyclic
dimers which correspond to two repeat units (or at least two potential
transesterification sites), the code must be modified before application
to systems containing two monomers. Care should also be taken to ensure
that only the reactions which are relevant to your particular system are
included. Due to the large number of parameters which are being
simulated, it is important that reasonable initial guesses are obtained
for all parameters.

The command line input order for both codes is: [exe_name] [input_name]
[output_name] [I0] [A0] [B0] [kAA] [kAB] [kBB] [kBA] [kAmA] [kAmB]
[kBmB] [kBmA] [kT] [accuracy_goal]

[input_name] and [output_name] are the names of the input and output
files, and are both read into the code as strings [I0], [A0], and [B0]
are the concentrations in mol/L or initiator, monomer A, and monomer B,
respectively. For the PLGA code, there is no difference between A and B,
both represent cyclic dimers. For the PLCL code, A is a cyclic dimer
(lactide in our case) and B is a monomer (caprolactone in our case) Each
of the k’s are initial guesses for the rate constants [accuracy goal]
sets the tolerance for SSR in the PLGA code. The code will continue to
run until the target SSR is reached (or indefinitely, if the goal is too
strict). We found that it was most effective to start with a large SSR,
and use the results from an initial fit with high SSR as inputs to then
fit again with lower SSR tolerance. We also found that eventually, the
simplex becomes too small to continue moving, and gets stuck at one
parameter set. It is of course possible to change the stop condition to
be when the simplex becomes sufficiently small, to avoid wasting time on
this issue. This change was implemented in the PLCL code, so in that
case [accuracy_goal] represents the average percent difference between
the worst point in the simplex and each of the other points. SSR is
still reported in the output file.

Running PLGA example: Command line input which can be used to generate
an output file similar to PLGA75-25.out: [exe_name] inputPLGA7525.in
[output_name].out 0.0405833 1.9737 5.9211 20.5 4 4 10 0.5 0.5 0.5 0.5
0.001 0.014

Since the model is stochastic, we can almost guarantee that the numbers
in your output file will not be exactly the same as ours! However, they
should be close. We ran the code 10 times with the command line input
shown above, and found the following standard deviations for each of the
rate constants.

kAA: 0.11, kAB: 0.33, kBB: 0.30, kBA: 0.21, kAmA: 0.19, kAmB: 0.40, 
kBmB: 0.16, KBmA: 0.20, kT: 0.0037

You can probably expect that the rate constants in your output file are
within a standard deviation or two of the ones in our output file.

Frequent troubleshooting steps: We have found that the acceptable
convention for defining the push_back function varies between compilers.
A namespace error upon compiling usually comes from this function (line
91 in fit-stochastic-reversible_di_v1.cpp and line 92 in
fit-stochastic-reversible_di_mo_v1.cpp)
