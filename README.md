# MITAR
Repository containing dynamical models for the MITAR project. 

R-scripts to compare pair-formation models of conjugation with bulk-conjugation models,
used in Alderliesten JB, Zwart MP, de Visser JAGM, Stegeman A, Fischer EAJ. 2022.
Second compartment widens plasmid invasion conditions: two-compartment pair-formation
model of conjugation in the gut. Journal of Theoretical Biology 533:110937.
* pairformation.R 
* pairformation2comp.R

The first file is used to model conjugation in a single compartment, the second file
to model conjugation in two compartments, with migration between them.
Select one parameter set from the section 'Parameter values' at a time,
run the section 'Run simulations', and run the appropriate output section.
Results of the simulations are saved as Comma Separated Values (.csv) files.
Graphs are created and, if the variabe 'saveplots' in the section 'Plotting and
simulation options' is set to TRUE, saved as .png-files.
If the variable 'plotdataapproxbulk' is set to TRUE, the data used to
approximate the bulk-conjugation rates is plotted.

R-scripts and matlab-script for multispecies models of conjugation:
Alderliesten et al. In prep.
* multispecies.R
* multispeciespinnewspecies.R
* multispecies_bifur.R
* multispecies_distr.R
* multispecies_analytic.m

These scripts are still under development and consider different scenarios
regarding how the plasmid-bearing population is introduced: in the
least-abundant or most-abundant species that is already present at the
plasmid-free equilibrium (multispecies.R), or in a new species that was not
yet present at the plasmid-free equilibrium and grows slower than, equally fast
as, or faster than the already-present species (multispeciespinnewspecies.R).
Scripts multispecies_bifur.R and multispecies_distr.R are used to create some
of the figures. The matlab-script multispecies_analytic.m is used to obtain
analytical expressions for the eigenvalues and derive (in)stability criteria
in the two-species case.
