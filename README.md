# MITAR
Repository containing dynamical models for the
[MITAR project](https://mitar.sites.uu.nl/). 


## Workflow
Most files in this repository are [R-scripts](https://cran.r-project.org/).
[RStudio](https://posit.co/download/rstudio-desktop/) is a useful IDE for R.
Help and tools to install R, R packages, and RStudio are available from the
[InstallPkgs GitHub repository](https://github.com/JesseAlderliesten/InstallPkgs).

The R scripts should be in the same folder as the `R Project (.Rproj)` file. The
used R-packages are loaded in a section `Loading required packages` at the top
of the script. The main packages are
[`deSolve`](https://CRAN.R-project.org/package=deSolve) (integrating ODEs),
[`rootSolve`](https://CRAN.R-project.org/package=rootSolve) (assess equilibria),
[`TruncatedNormal`](https://CRAN.R-project.org/package=TruncatedNormal) (random
number generation from truncated distributions),
[`dplyr`](https://CRAN.R-project.org/package=dplyr) for data handling, and
[`ggplot2`](https://CRAN.R-project.org/package=ggplot2),
[`scales`](https://CRAN.R-project.org/package=scales) (to display data and results).


## File overview

R-scripts to compare pair-formation models of conjugation with bulk-conjugation
models, used in Alderliesten JB, Zwart MP, de Visser JAGM, Stegeman A, Fischer
EAJ. 2022. Second compartment widens plasmid invasion conditions:
two-compartment pair-formation model of conjugation in the gut. Journal of
Theoretical Biology 533:110937.
* pairformation.R 
* pairformation2comp.R

The first script is used to model conjugation in a single compartment, the
second script to model conjugation in two compartments with migration between
them. Select one parameter set from the section 'Parameter values' at a time,
run the section 'Run simulations', and run the appropriate output section.
Results of the simulations are saved as Comma Separated Values (.csv) files.
Graphs are created and, if the variabe 'saveplots' in the section 'Plotting and
simulation options' is set to TRUE, saved as .png-files. If the variable
'plotdataapproxbulk' is set to TRUE, the data used to approximate the
bulk-conjugation rates is plotted.

R-scripts and MATLAB-script for multispecies models of conjugation:
Alderliesten JB, Zwart MP, de Visser JAGM, Stegeman A, Fischer EAJ. In prep.
The effects of ecological interactions and distinct conjugation rates on the
invasion of a conjugative plasmid in bacterial communities.

These scripts use a generalised Lotka-Volterra model elaborated with
plasmid-bearing populations and conjugation to investigate how the richness,
evenness, and diversity of a plasmid-free microbiome influence the ability of
plasmid-bearing bacteria to invade it (results for plasmid-free bacteria as
invaders are also included in the results of the simulations). Throughout, the
interspecies conjugation rates involving the invading species (i.e., the
initially plasmid-bearing species) are either equal to the other interspecies
conjugation rates or are thousandfold lower.

* multispecies.R: simulate invasion through a species that is already present at
the plasmid-free equilibrium. Interaction coefficients are drawn from
distributions, species abundandances are specified by species abundance models,
and then growth rates are calculated as r* = -AB* to obtain an equilibrium.
* multispecies_altgrow.R: similar to multispecies.R, but here both interaction
coefficients and growth rates are drawn from distributions, and then species
abundances are calculated as B* = - A^-1 r to obtain an equilibrium.
* multispecies_analytic.m: MATLAB-script to obtain an analytical expressions for
the eigenvalues and derive (in)stability criteria in the two-species case.
* multispecies_bifur.R: R-script to create bifurcation-like plots showing how
the intraspecies conjugation rate of the initially plasmid-bearing species and
fitness costs of bearing a plasmid influence the epidemiological stability of
the plasmid-free equilibrium (Figure 2 and Figure S6).
* multispecies_distr.R: R-script to create examples of probability densities for
the distributions of ecological interaction coefficients (Figure S1).
* multispecies_funcs.R: R-script containing the functions used by the other
R-scripts.
* multispecies_pinnew.R: similar to multispecies.R, but here invasion through a
species that was not yet present at the plasmid-free equilibrium is simulated.
* multispecies_pinnew_bifur.R: similar to multispecies_bifur.R, but here
invasion through a species that was not yet present at the plasmid-free
equilibrium is simulated.

To run the simulations, choose one of the R-scripts and run it from the top to
load the required packages, select the basis parameter set from the section
'Settings and defining parameter space', run the simulations and plot the
results. Results of the simulations are saved as R-data (.RDS) files. Graphs are
created and, if the variabe 'saveplots' in the section 'Settings and defining
parameter space' is set to TRUE, saved as .png-files. If the variable
'saveplotconjovertime' is set to TRUE, plots showing the dynamics over time are
created and saved.

Changes to the settings allow to model further scenarios. For example, the
initial total abundance, fitness costs of bearing a plasmid, and interspecies
conjugation rate can be changed. In addition, the default scenario where
plasmid-bearing bacteria replace bacteria of the most-abundant species can be
changed to let them replace bacteria of the least-abundant species by changing
variable 'PReplMostAbun' from TRUE (the default) to FALSE. For scripts that
simulate invasion through a new species, the new species that introduces the
plasmid grows slower than (1), equally fast as (2, default), or faster than (3)
the already-present species, depending on the value of variable
'newgrowthratecode'.
