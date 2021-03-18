# MITAR
Repository containing dynamical models for the MITAR project

Files used to model pair formation are:
* pairformation.R
* pairformation2comp.R
These scripts are used to run simulations with the pair-formation models.
Select one parameter set from the section 'Parameter values' at a time,
run the section 'Run simulations', and run the appropriate output section.
Results of the simulations are saved as Comma Separated Values (.csv) files.
Graphs are created and, if the variabe 'saveplots' in the section 'Plotting and
simulation options' is set to TRUE, saved as .png-files.
If the variable 'plotdataapproxbulk' is set to TRUE, the data used to
approximate the bulk-conjugation rates is plotted.


Files to create a multispecies model of conjugation:
* multispecies.R
This model is still under development.

