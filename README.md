# MITAR
Repository containing dynamical models for the MITAR project

R-scripts to compare pair-formation models of conjugation with bulk-conjugation models:
* pairformation.R 
* pairformation2comp.R (m)

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
* multispeciespinleastabun.R
* multispeciespinnewspecies.R

These scripts are for three different scenarios regarding how the plasmid-bearing
population is introduced: in the most-abundant species, in the least abundant species,
or in a species that was not yet present model. These models are is still under development.

R-script to analyze microbiome data of the BEWARE project:
* bewaremicrobio.R (still under development)
