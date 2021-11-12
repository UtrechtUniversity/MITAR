# MITAR
Repository containing dynamical models for the MITAR project

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

R-scripts for multispecies models of conjugation:
* multispecies.R
* multispeciespinnewspecies.R

These scripts consider three different scenarios regarding how the plasmid-bearing
population is introduced: in the most-abundant species, in the least-abundant species,
or in a species that was not yet present. These models are still under development.
