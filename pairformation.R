#### Pair-formation model of conjugation ####

#### Version history ####
# Before 9 January 2020 the equation for dR contained - kp*R*R instead of - 2*kp*R*R
# On 10 february 2020 I changed the equation of the nutrients to include washout of nutrients
# Later I changed kp, kn, gd, and gt to log10kp, log10kn, log10gd, and log10gt
# In version 9 (13 March 2020) I removed the biomass equations from the models. This has the
# advantage that the integration will be slightly faster since one equation less has to be
# integrated. More importantly, removing the biomass equation is needed to enable
# stability-analysis based on the eigenvalues of the jacobimatrix in the plasmid-free equilibrium
# On 13 March I also removed the plotting commands from the script and changed the structure of
# the output dataframe, to prevent repetition of the biomass and Dinit = plasmid-bearing bacteria at the initial state.
# For some versions there is a script available using rootSolve to run fast to equilibrium,
# e.g. PairFormationSimulation2rootSolve8, but no plot over time can be made and it does not always work,
# and the output structure is different from deSolve, leading to different selection styles for storing data
# in the dataframe.
# In version 13 I corrected the calculation of the eigenvalues to use the plasmid-free equilibrium EqFull, not
# the plasmid-free equilibrium with Dinit donors added to it (state). I also changed the calculation of the
# biomass at the initial state in the dataframe to a cleaner implementation.
# Version 13 is a copy of version 12, where I copied the script part by part to get the correct outlining
# of all the brackets, and I changed some comments and the order of some parts to get a clearer and cleaner
# script.
# In version 13 (?) added dots in the graphs to indicate the points where equilibrium was not reached.
# Besides, I splitted cSet and gSet in separate sets for D and T, and added Dinit as a set and adjusted the
# for-loops accordingly.
# In version 14 I added automated labelling of the facets, and changed the manual labelling of the x- and
# y-axis to automated labelling.
# Versions 15 - 17 exist in a version with complete pair-formation, and a version ('ChangedContact') where
# Zhong's model is followed regarding pair-formation, such that only Mdr and Mrt pairs are formed, and Mdt
# and Mtt pairs arise from conjugation
# In version 18 I changed the linestyles and linecolors to highlight the plasmid-free population in the same
# way as in the earlier versions of the 'complete' model.
# In version 19 I extended the root-function and added an eventfunction to set populations of bacteria that are
# smaller than a threshold to 0. I also added handling of ratios that are invalid when checking for coexistence.
# In version 21 I changed the order of equations to have nutrients as the first equation so nutrients 
# can be disregarded in the rootfunction by using state[-1] without having to call length(state) or
# length(state) + 1. Similarly I changed the order of the roots in the rootfunction, to have terminalroot = 1
# instead of length(state) + 1. I also changed the rootfunction by selecting for only the bacterial populations
# to trigger a root. After these changes I also changed mylty, mycol, and initial states to match the changed order.
# NOTE: this implementation of events was WRONG, because the nutrients were excluded from selection so the
# wrong states were put to 0. This issue was solved in later versions of version 21.
# In version 23 the calculation of the bulk-rates according to Zhong has been added,
# and gdbulk and gtbulk have been added to the parameters stored in the matrix.
# The root- and eventfunctions are now also incorporated into the bulk-conjugation model.
# The the bulk-conjugation model is ran using Zhong's bulk-conjugation rates, and
# the output is stored in the matrix as well. I have added the column ComplexEigVal
# to the matrix, which indicates if complex eigenvalues were present, and now only
# the real parts of the eigenvalues are stored. Previously all values in the matrix
# would have a complex +0i part added if a value with a complex part was stored,
# and the imaginary parts hamper other functions.
# In version 24 I changed to writing csv-files instead of xlsx files, changed the
# script such that, even if invasion with the pair-formation model is not possible,
# gdbulk and gtbulk are calculated, and stability-analysis of the plasmid-free
# equilibrium in the bulk-model (which I added) is performed.
# I also lowered extinctionthreshold to prevent the density of newly formed
# populations from being below the extinctionthreshold.
# In version 25 I changed NI from 1 to 10 (which has previously been used),
# to decrease simulation time. The plots were also updated, various new plots
# were added, and other colorscales were used.
# In version 26 I removed the loop for e, Mytmax, Mytstep, and smallchange.
# The values for kp and kn were changed to have integer values for log10(kp)
# and log10(kn) (they were seq(from = -10.5, to = -8.5, by = 0.5), but the graphs
# erroneously present it like they were -10, -9, ...) 
# The loops for kp and kn were moved to only change after the plasimd-free
# equilibrium has been calculated.  


#### To do ####
# Add reference section for references to help-files, powerpoints, publications etc.
# Try to write plotting function with vars to be plotted on x- and y-axis and in facets,
# difference between models or ratio between models as options.

# Add totalT bulk, totalR bulk etc at equilibrium to dataframe
# Run bulk-model for short time (same as short pair-model) and compare them
# to see if the bulk-conjugation parameter holds (automated analysis somehow)

# Remove more out of loop?
# The runs which reach tmax were finished very quickly when I ran them
# out of the loop (manually setting the variables to the right numbers).

# QUESTION: have the resources to be considered for stability analysis?
# Imran does exclude them from stability analysis?
# Since I am only interested in perturbations of the equilibrium by adding donors,
# can the stability analysis be simplified (by only looking at the eigenvectors
# (eigenvalue?) of the donor-equation?). See 'A Practical Guide to Ecological
# Modelling: Using R as a Simulation Platform', p. 229-230)

#### Introduction ####
# The idea for such a model was taken from equation (2) of the following publication:
# Zhong X, Krol JE, Top EM, Krone SM. 2010. Accounting for mating pair formation
# in plasmid population dynamics. Journal of Theoretical Biology 262:711-719.

# See also: Zhong X, Droesch J, Fox R, Top EM, Krone SM. 2012. On the meaning and
# estimation of plasmid transfer rates for surface-associated and well-mixed
# bacterial populations. Journal of Theoretical Biology 294:144-152.

# To read data from csv-file
# FileName <- "2020_mei_28_10_outputdeSolveChangedContact_Samengevoegd.csv"
# MyData <- read.csv(FileName, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE
# )
# MyData <- as.data.frame(MyData)



#### Model functions ####
# Model a population of plasmid-free recipients with nutrients. I do not use
# this model in the script, because instead I use the analytical solution for
# the plasmid-free equilibrium as found with Matlab.
ModelRecipNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w - e*bR*Nutr*R
    dR <- bR*Nutr*R - w*R
    return(list(c(dNutr, dR)))
  })
}

# The pair-formation conjugation model. Mdr pairs and Mrt pairs are formed with
# rate kp*R if recipients meet donors or transconjugants, respectively.
# Conjugation occurs in these Mdr and Mrt pairs, with intrinsic conjugation rate
# gd and gt, respectively. From conjugation, Mdt and Mtt pairs arise. The
# different pairs fall apart again with rate kn. In this structure of pair-formation 
# I follow Zhong's model. As addition to Zhong's model, I included costs, washout,
# and nutrients.
ModelPairsNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w - e*Nutr*((1 - cd)*bR*(D + Mdr + Mdt) + bR*(R + Mdr + Mrt) +
                                       (1 - ct)*bR*(Trans + Mdt + Mrt + 2*Mtt))
    dD <- (1 - cd)*bR*Nutr*(D + Mdr + Mdt) - 10^(log10kp)*D*R + 10^(log10kn)*(Mdr + Mdt) - w*D
    dR <- bR*Nutr*(R + Mdr + Mrt) - 10^(log10kp)*R*(D + Trans) + 10^(log10kn)*(Mdr + Mrt) - w*R
    dTrans <- (1 - ct)*bR*Nutr*(Trans + Mdt + Mrt + 2*Mtt) - 10^(log10kp)*R*Trans + 10^(log10kn)*(Mdt + Mrt + 2*Mtt) -
      w*Trans
    dMdr <- 10^(log10kp)*D*R - 10^(log10kn)*Mdr - 10^(log10gd)*Mdr - w*Mdr
    dMdt <- 10^(log10gd)*Mdr - 10^(log10kn)*Mdt - w*Mdt
    dMrt <- 10^(log10kp)*R*Trans - 10^(log10kn)*Mrt - 10^(log10gt)*Mrt - w*Mrt
    dMtt <- 10^(log10gt)*Mrt - 10^(log10kn)*Mtt - w*Mtt
    return(list(c(dNutr, dD, dR, dTrans, dMdr, dMdt, dMrt, dMtt)))
  })
}

# Model to approximate the bulk-conjugation rate of the transconjugant.
# Nutrients, growth, washout, and donors are not included in this model.
ModelEstConjBulkTrans <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR <- - 10^(log10kp)*R*Trans + 10^(log10kn)*Mrt
    dTrans <- - 10^(log10kp)*R*Trans + 10^(log10kn)*(Mrt + 2*Mtt)
    dMrt <- 10^(log10kp)*R*Trans - 10^(log10kn)*Mrt - 10^(log10gt)*Mrt
    dMtt <- 10^(log10gt)*Mrt - 10^(log10kn)*Mtt
    return(list(c(dR, dTrans, dMrt, dMtt)))
  })
}

# The bulk-conjugation model
ModelBulkNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w - e*Nutr*((1 - cd)*bR*D + bR*R + (1 - ct)*bR*Trans)
    dD <- (1 - cd)*bR*Nutr*D - w*D
    dR <- bR*Nutr*R - gdbulk*D*R - gtbulk*Trans*R - w*R
    dTrans <- (1 - ct)*bR*Nutr*Trans + gdbulk*D*R + gtbulk*Trans*R - w*Trans
    return(list(c(dNutr, dD, dR, dTrans)))
  })
}
# The plasmid-free equilibrium is R* = (NI - (w/bR))/e, Nutr* = w/bR, D* = Trans* = 0
# with condition w < NI*bR (otherwise it becomes a cell-free equilibrium ?) and
# the stability of this equilibrium is determined by the following eigenvalue
# - ct*w - (gt*(w - NI*bR))/(bR*e), leading to gt > (bR*ct*e*w)/(NI*bR - w) as
# invasion criterion. THIS SUGGESTS THAT DONOR CHARACTERISTICS DO NOT PLAY A ROLE
# in determining invasion

# Define root-function and event-function
# If the root-argument is equal to 0, an event (as defined by the event-function)
# is triggered. See ?events and ?lsodar (both in the deSolve package) for
# background information and examples. 
# The first root becomes 0 if the sum of absolute rates of change is equal to the
# threshold smallchange, i.e., when equilibrium is nearly reached. This root is
# set as the terminal root, in order to terminate the integration if this occurs.
# The other roots are used to determine if any of the state variables is equal
# to the extinction threshold. If this occurs, the event function is called,
# which sets the states that are equal to threshold to zero, and the integration
# continues.
# NOTE that the NUTRIENTS are also in the root- and event-functions, so currently
# they would also be set to 0 if they fall below the threshold, which is undesirable.
# The nutrients can probably removed from the rootfunction by using state[-1] -
# extinctionthreshold, but attempts to remove them from events failed for events
# implemented as state2[state2[-1] < extinctionthreshold] <- 0 or state2[state2[-1]
# < extinctionthreshold] <- 0 because the index shifts if nutrients are below the
# threshold. Using state2[(state2[-1] < extinctionthreshold) + 1] <- 0 might work.
# Note that it is still possible to have bacterial densities lower than the
# threshold, if these populations are created from other populations at each
# timestep.
rootfun <- function(t, state, parmspair) {
  c(sum(abs(unlist(ModelPairsNutr(t, state, parmspair)))) - smallchange, state - extinctionthreshold)
}

rootfunBulk <- function(t, stateBulk, parmsBulk) {
  c(sum(abs(unlist(ModelBulkNutr(t, stateBulk, parmsBulk)))) - smallchange, stateBulk - extinctionthreshold)
}

eventfun <- function(t, state, parmspair) {
  state[state < extinctionthreshold] <- 0
  return(state)
}

eventfunBulk <- function(t, stateBulk, parmsBulk) {
  stateBulk[stateBulk < extinctionthreshold] <- 0
  return(stateBulk)
}

#### Explanation of parameter values ####
# Growthrates (1/h) for recipient: bR (growth rates for donor and transconjugant
# are not explicitly defined, but calculated as (1 - cd)*bR and (1 - ct)*bR)
# Conjugation rates (1/h) from donor and transconjugant: gd and gt
# Mating pair attachment rate (mL * cell^-1 * h^-1): kp (k+ in Zhong's notation)
# Mating pair detachment rate (1/h): kn (k- in Zhong's notation)
# Nutrient concentration in the inflowing liquid (microgram * mL^-1): NI
# Resource conversion rate (microgram per cell division): e
# Washout rate (1/h): w

### Parameter values
# Stable co-existence of recipients, transconjugants, and Mrt and Mtt pairs
# time = 816943.9, Nutr = 0.02427729 D = 0, R = 3829218, Trans = 6142815, Mdr = Mdt = 0, Mrt = 324.6586, Mtt = 1520.435
# bRSet <- c(1.7)
# NISet <- c(10)
# log10kpSet <- c(-9.6)
# log10knSet <- c(0.5)
# eSet <- c(1E-6)
# wSet <- c(0.04)
# cdSet <- c(0.05)
# ctSet <- c(0.05)
# log10gdSet <- c(1.176)
# log10gtSet <- c(1.176)

# Single run, invasion not possible
DinitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
wSet <- c(0.04)
log10kpSet <- c(-10)
log10knSet <-  c(0.3)
cdSet <- c(0.05)
ctSet <- c(0.05)
log10gdSet <- c(1.176)
log10gtSet <- c(1.176)
NutrConv <- c(1E-6)


# Single run, invasion possible
DinitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
wSet <- c(0.04)
log10kpSet <- c(-6)
log10knSet <-  c(0.3)
cdSet <- c(0.05)
ctSet <- c(0.05)
log10gdSet <- c(1.176)
log10gtSet <- c(1.176)
NutrConv <- c(1E-6)

# Vary kp
DinitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
wSet <- c(0.04)
log10kpSet <- seq(from = -11, to = -5, by = 0.5)
log10knSet <- c(0.3)
cdSet <- c(0.05)
ctSet <- c(0.05)
log10gdSet <- c(1.176)
log10gtSet <- c(1.176)
NutrConv <- c(1E-6)

# # Vary kp and kn
# DinitSet <- c(1E3)
# bRSet <- c(1.7)
# NISet <- c(10)
# wSet <- c(0.04)
# log10kpSet <- seq(from = -11, to = -5, by = 0.5)
# log10knSet <- seq(from = -1, to = 3, by = 0.5)
# cdSet <- c(0.05)
# ctSet <- c(0.05)
# log10gdSet <- c(1.176)
# log10gtSet <- c(1.176)
# NutrConv <- c(1E-6)

# # Vary kp, kn, cd, and ct
# DinitSet <- c(1E3)
# bRSet <- c(1.7)
# NISet <- c(10)
# wSet <- c(0.04)
# log10kpSet <- seq(from = -11, to = -5, by = 1)
# log10knSet <- seq(from = -1, to = 3, by = 1)
# cdSet <- c(0.01)
# ctSet <- c(0.01, 0.05)
# log10gdSet <- c(1.176)
# log10gtSet <- c(1.176)
# NutrConv <- c(1E-6)

# # Vary kp, kn, gd, gt, cd, and ct
DinitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
wSet <- c(0.04)
log10kpSet <- seq(from = -11, to = -5, by = 0.1)
log10knSet <- seq(from = -1, to = 3, by = 0.1)
cdSet <- c(0.01, 0.025, 0.05)
ctSet <- c(0.01, 0.025, 0.05)
log10gdSet <- c(1, 1.176)
log10gtSet <- c(1, 1.176)
NutrConv <- c(1E-6)

# # Extensive dataset
# DinitSet <- c(1E3)
# bRSet <- c(0.4, 1.7)
# NISet <- c(0.1, 1, 10)
# wSet <- c(0.04)
# log10kpSet <- seq(from = -11, to = -8, by = 0.5)
# log10knSet <- seq(from = -1, to = 3, by = 0.5)
# cdSet <- c(0.01, 0.025, 0.05)
# ctSet <- c(0.01, 0.025, 0.05)
# log10gdSet <- c(1, 1.176)
# log10gtSet <- c(1, 1.176)
# NutrConv <- c(1E-6)


DinitSet <- c(1E3)
log10kpSet <- seq(from = -11, to = -5, by = 1)
log10knSet <- seq(from = -1, to = 3, by = 1)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
log10gdSet <- c(1, 1.176)
log10gtSet <- c(1, 1.176)

#### Create matrix to store data ####
Mydf <- expand.grid(Dinit = DinitSet, bR = bRSet, NI = NISet, w = wSet, log10kpSet = log10kpSet,
                    log10knSet = log10knSet, cdSet = cdSet, ctSet = ctSet,
                    log10gdSet = log10gdSet, log10gtSet = log10gtSet, KEEP.OUT.ATTRS = FALSE)


TotalIterations <- length(bRSet)*length(NISet)*length(log10kpSet)*
  length(log10knSet)*length(wSet)*length(DinitSet)*length(cdSet)*length(ctSet)*
  length(log10gdSet)*length(log10gtSet)
print(TotalIterations)
CurrentIteration <- 0
MyData <- matrix(data = NA, nrow = TotalIterations, ncol = 61, byrow = TRUE) # To store output data
DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")

# Instead of performing everything inside the loop, one could also use the loop
# to create paramtervalues, and then for each row in the matrix perform the
# simulations, or use a list directly, as in the following example
# outlist <- list()
# plist <- cbind(r = runif(10, min = 0.1, max = 5),
#                K = runif(10, min = 8, max = 15))
# for (i in 1:nrow(plist))
#   outlist[[i]] <- ode(y = c(y = 2), times, derivs, parms = plist[i,])
# plot(out, outlist)


# Times for which output of the simulation is wanted.
# Note that I don't use timestep untill t = 100, and note furthermore that
# the used ode-solvers are variable-step methods, so the times in times
# are NOT the only times at which integration is performed. See
# ?diagnostics.deSolve() and ?lsodar() for details.
if(runsimulation == 1) {
  times <- c(0:100, seq(from = 100 + Mytstep, to = Mytmax, by = Mytstep)) 
}

################################################################################

library(deSolve) # To solve differential equations
library(rootSolve) # To get jacobian matrix for stability-analyses
library(tidyr) # for 'expand.grid()' with dataframe as input
library(dplyr)
library(ggplot2) # For plotting data
library(RColorBrewer) # For better color schemes

#### Plotting options ####
mylty <- c(lty = c(3, 1, 2, 1, 1, 1, 1, 1))
# mycol <- c("black", "purple", "hotpink", "red", "yellow", "green1", "blue", "cyan")
mycol <- c("black", brewer.pal(7, "Set1"))
MyColorBrew <- brewer.pal(11, "Spectral") # see display.brewer.all()
MyColorBrew2 <- brewer.pal(9, "YlOrRd")
timesParmsEst <- seq(from = 0, to = 3, by = 0.1)
myylim <- c(1E-4, 1E7) # Defining the limits for the y-axis
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
runsimulation <- 0 # if runsimulation == 0, no simulation is run, and only the
# parametervalues, plasmid-free equilibrium, and the eigenvalues are stored
plotoutput <- 0
extinctionthreshold <- 1E-10 # Population size is set to 0 if it is below the extinctionthreshold
verbose <- 0 # if verbose == 1, diagnositics on the simulations are printed and roots are indicated in the graphs
smallchange <- c(1E-5)
Mytmax <- c(1E5)
Mytstep <- c(10)
# Model to approximate the bulk-conjugation rate of the donor.
# Nutrients, growth, washout, conjugation from transconjugants, and Mtt-pairs
# are not included in this model.
ModelEstConjBulkDonor <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dD <- - 10^(log10kp)*D*R + 10^(log10kn)*(Mdr + Mdt)
    dR <- - 10^(log10kp)*R*(D + Trans) + 10^(log10kn)*(Mdr + Mrt)
    dTrans <- - 10^(log10kp)*R*Trans + 10^(log10kn)*(Mdt + Mrt)
    dMdr <- 10^(log10kp)*D*R - 10^(log10kn)*Mdr - 10^(log10gd)*Mdr
    dMdt <- 10^(log10gd)*Mdr - 10^(log10kn)*Mdt
    dMrt <- 10^(log10kp)*R*Trans - 10^(log10kn)*Mrt
    return(list(c(dD, dR, dTrans, dMdr, dMdt, dMrt)))
  })
}

## Small parameterset for tests
bRSet <- c(1.7)
NISet <- c(10, 100)
wSet <- c(0.04)
NutrConv <- c(1E-6)
log10kpSet <- c(-10, -6)
log10knSet <- c(0.3)
cdSet <- c(0.05)
ctSet <- c(0.01, 0.05)
log10gdSet <- c(1.176)
log10gtSet <- c(1.176)
DinitSet <- 1000

## Large dataset for tests
DinitSet <- c(1E3)
bRSet <- c(0.6, 1.7)
NISet <- c(10, 100)
wSet <- c(0.04)
NutrConv <- c(1e-6)
log10kpSet <- seq(from = -11, to = -5, by = 1)
log10knSet <- seq(from = -1, to = 3, by = 1)
log10gdSet <- c(1, 1.176)
log10gtSet <- c(1, 1.176)


## Calculate plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(bR = bRSet, NI = NISet, NutrConv = NutrConv, w = wSet)
if(any(MyData <= 0)) warning("All parameters should have positive values.")

# Calculate the plasmid-free equilibrium (R*, Nutr*)

calceqplasmidfree <- function(MyData) {
  with(as.list(MyData), {
    REq <- ((NI - w / bR)) / NutrConv
    NutrEq <- w / bR
    Eq <- c(NutrEq = NutrEq, REq = REq)
    return(Eq)
  })
}

dfeqplasmidfree <- apply(X = MyData, MARGIN = 1, FUN = calceqplasmidfree)
dfeqplasmidfree <- t(dfeqplasmidfree)     # ToDo: change the function to get the transpose as output

# Check if plasmid-free equilibrium is positive, warn if not.
if(any(dfeqplasmidfree[, "REq"] <= 0)) warning("The number of recipients at equilibrium is not positive!
Increase the nutrient concentration in the inflowing liquid by changing NI?")

MyData <- cbind(MyData, dfeqplasmidfree)

# Note on indexing of tibbles: use MyData[[1, "REq"]] to return a vector (MyData[1, "REq"] returns a LIST)
# See is.vector(state) and is.numeric(state)
# Using MyData$Dinit[i] in a loop does work with tibbles as well

# Add combinations with the parameters needed to approximate gdbulk and gtbulk to MyData
MyData <- expand_grid(MyData, log10kp = log10kpSet, log10kn = log10knSet,
                       log10gd = log10gdSet, log10gt = log10gtSet, Dinit = DinitSet)

# ToDo: try to find method to obtain the order in which I specify the state from EstConjBulkDonor
# or from ModelEstConjBulkDonor, to prevent hardcoding names on the returned object.

# Run simulation with adjusted pair-formation models for a short timespan and
# calculate approximations of gdbulk and gtbulk from the output at t = 3 hours,
# following Zhong's approach for the calculations. NOTE: I used timesteps of 0.1
# instead of 1, because 3 timesteps is too few to get stable estimates. I did
# not use root- or eventfunctions here because they lead to unstable behaviour
# (in the early timesteps if invasion is possible, or over the whole simulation
# if invasion is not possible), leading to unstable estimates of gdbulk and
# gtbulk.
EstConjBulkDonor <- function(MyData) {
  state <- c(D = MyData[["Dinit"]], R = MyData[["REq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
  parms <- MyData
  DataEstConjBulkDonor <- tail(ode(t = timesParmsEst, y = state,
                                   func = ModelEstConjBulkDonor, parms = parms), 1)
  #print("DataEstConjBulkDonor=")
  #print(DataEstConjBulkDonor)
  return(DataEstConjBulkDonor)
}

DataEstConjBulkDonor <- apply(X = MyData, MARGIN = 1, FUN = EstConjBulkDonor)
DataEstConjBulkDonor <- t(DataEstConjBulkDonor)
colnames(DataEstConjBulkDonor) <- c("time", "D", "R", "Trans", "Mdr", "Mdt", "Mrt")

TotalDEstConjBulkDonor <- DataEstConjBulkDonor[, "D"] + DataEstConjBulkDonor[, "Mdr"] + DataEstConjBulkDonor[, "Mdt"]
TotalREstConjBulkDonor <- DataEstConjBulkDonor[, "R"] + DataEstConjBulkDonor[, "Mdr"] + DataEstConjBulkDonor[, "Mrt"]
gdbulk <- (10^MyData[, "log10gd"]) * DataEstConjBulkDonor[, "Mdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor)
gdbulk <- unname(gdbulk)
MyData <- cbind(MyData, gdbulk = gdbulk)
MyData

# This works and results of gdbulk are identical to those of earlier versions of the script containing the nested for-loops.
write.csv(MyData, file = "Calculategdbulkapply.csv", quote = FALSE, row.names = FALSE)





DataEstConjBulkDonor <- tail(ode(t = timesParmsEst, y = stateDonor,
                                 func = ModelEstConjBulkDonor, parms = parmsEstConjBulk), 1)
TotalDEstConjBulkDonor <- DataEstConjBulkDonor[, "D"] + DataEstConjBulkDonor[, "Mdr"] + DataEstConjBulkDonor[, "Mdt"]
TotalREstConjBulkDonor <- DataEstConjBulkDonor[, "R"] + DataEstConjBulkDonor[, "Mdr"] + DataEstConjBulkDonor[, "Mrt"]
gdbulk <- (10^log10gdValue) * DataEstConjBulkDonor[, "Mdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor)

stateTrans <- c(R = RAna1, Trans = Dinit, Mrt = 0, Mtt = 0)
DataEstConjBulkTrans <- tail(ode(t = timesParmsEst, y = stateTrans,
                                 func = ModelEstConjBulkTrans, parms = parmsEstConjBulk), 1)
TotalTEstConjBulkTrans <- DataEstConjBulkTrans[, "Trans"] + DataEstConjBulkTrans[, "Mrt"] + 2*DataEstConjBulkTrans[, "Mtt"]
TotalREstConjBulkTrans <- DataEstConjBulkTrans[, "R"] + DataEstConjBulkTrans[, "Mrt"]
gtbulk <- (10^log10gtValue) * DataEstConjBulkTrans[, "Mrt"] / (TotalREstConjBulkTrans * TotalTEstConjBulkTrans)


                
                
                # # Store equilibria for determination of stability, add donors 
                # # to make state the starting point for simulations.
                # 
                # EqFull <- c(Nutr = NutrAna, D = 0, R = RAna1, Trans = 0,
                #             Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0) 
                # state <- EqFull
                # state["D"] <- Dinit
                # EqFullBulk <- c(Nutr = NutrAna, D = 0, R = RAna1, Trans = 0)
                # stateBulk <- EqFullBulk
                # stateBulk["D"] <- Dinit
                
                
                
                
                for(cdValue in cdSet) {
                  for(ctValue in ctSet) {

                    CurrentIteration <- CurrentIteration + 1
                    if(round(CurrentIteration / 500) == CurrentIteration / 500) {
                      print(paste0("Current iteration = ", CurrentIteration,
                                   ", time = ",
                                   format(Sys.time(), format = "%H:%M:%S")))
                    }
                    
                    parmspair <- c(bR = bRValue, NI = NIValue, e = NutrConv, w = wValue,
                                   log10kp = log10kpValue, log10kn = log10knValue, 
                                   cd = cdValue, ct = ctValue, log10gd = log10gdValue,
                                   log10gt = log10gtValue)
                    
                    parmsBulk <- c(parmspair, gdbulk = gdbulk, gtbulk = gtbulk)
                    
                    # Numerically estimate the Jacobian matrix of the plasmid-free
                    # equilibrium of the pair-formation model, then calculate
                    # (or approximate?) the eigenvalues of this matrix.
                    EigValEq <- eigen(x = jacobian.full(y = EqFull, func = ModelPairsNutr, parms = parmspair),
                                      symmetric = FALSE, only.values = TRUE)$values
                    ComplexEigVal <- is.complex(EigValEq) 
                    EigValEq <- Re(EigValEq) # selecting real part of the eigenvalues
                    DomEigVal <- max(EigValEq) # selecting the maximum real part of the eigenvalues
                    SignDomEigVal <- sign(DomEigVal)
                    SignEigValEqual <- identical(rep(SignDomEigVal, length(EigValEq)), sign(Re(EigValEq)))
                    
                    # Determine eigenvalues of the jacobian matrix of the plasmid-free equilibrium
                    # of the bulk-conjugation model. If only.values = FALSE, the eigenvectors are
                    # stored in $vectors. See ?jacobian.full() for an example.
                    EigValEqBulk <- eigen(x = jacobian.full(y = EqFullBulk, func = ModelBulkNutr,
                                                            parms = parmsBulk),
                                          symmetric = FALSE, only.values = TRUE)$values
                    ComplexEigValBulk <- is.complex(EigValEqBulk) 
                    EigValEqBulk <- Re(EigValEqBulk) # selecting real part of the eigenvalues
                    DomEigValBulk <- max(EigValEqBulk) # selecting the maximum real part of the eigenvalues
                    SignDomEigValBulk <- sign(DomEigValBulk)
                    SignEigValEqualBulk <- identical(rep(sign(DomEigValBulk), length(EigValEqBulk)),
                                                     sign(Re(EigValEqBulk)))
                    
                    if(SignDomEigVal == -1) {
                      # No invasion possible in the pair-formation model
                      EqAfterInvDonor <- c(time = 0, EqFull)
                    } else {
                      if(runsimulation == 1) {
                        
                        # Invasion is possible, run simulation to see how
                        # many plasmid-bearing bacteria are present at equilibrium
                        
                        # The initial state for phase 2 is the plasmid-free
                        # equilibrium (R*, Nutr*) with the addition of
                        # Dinit donor bacteria per mL
                        out2 <- ode(t = times, y = state, func = ModelPairsNutr,
                                    parms = parmspair, rootfun = rootfun,
                                    events = list(func = eventfun, root = TRUE, terminalroot = 1))
                        if(verbose == TRUE) {
                          print(diagnostics(out2))
                          print(attributes(out2))
                        }
                        EqAfterInvDonor <- tail(out2, 1)[, ]
                        
                        if(plotoutput == 1) {
                          
                          Ratio <- EqAfterInvDonor["Trans"]/(EqAfterInvDonor["Trans"] + EqAfterInvDonor["R"])
                          if(is.na(Ratio)==FALSE) {
                            if(Ratio > 0.001 & Ratio < 0.999) {
                              Coexistence <- 1
                            } else {
                              Coexistence <- 0
                            }
                          } else {
                            Coexistence <- 0
                          }
                          # if(Coexistence == 1) {
                          maintitle <- c("Pair-formation model")
                          subtitlepair <- paste0("bR=", bRValue,
                                                 " NI=", NIValue,
                                                 " log10kp=", log10kpValue,
                                                 " log10kn=", log10knValue,
                                                 " e=", NutrConv,
                                                 " w=", wValue,
                                                 " cd=", cdValue,
                                                 " ct=", ctValue,
                                                 " log10gd=", log10gdValue,
                                                 " log10gt=", log10gtValue)
                          # png(filename = paste0(DateTimeStamp, "longrun", CurrentIteration + 1, ".png"))
                          if(verbose == TRUE) {
                            matplot.deSolve(out2, main = maintitle,
                                            sub = subtitlepair, ylim = myylim,
                                            xlim = c(0, tail(attributes(out2)$troot, 2)[1]),
                                            log = if(yaxislog == 1) {"y"},
                                            col = mycol, lty = mylty, lwd = 2,
                                            legend = list(x = "topright"))
                            abline(h = extinctionthreshold)
                            abline(v =  attributes(out2)$troot)
                            grid()
                          } else {
                            matplot.deSolve(out2, main = maintitle, sub = subtitlepair, ylim = myylim,
                                            log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                                            legend = list(x = "bottomright"))
                            grid()                                      
                          }
                          
                          
                          # dev.off()
                          # }
                        } # End of preparing and plotting
                      } # End of simulation part
                    }
                    
                    if(SignDomEigValBulk == -1) {
                      # No invasion possible in the bulk-conjugation model
                      EqAfterInvDonorBulk <- c(time = 0, EqFullBulk)
                    } else {
                      if(runsimulation == 1) {
                        # Run bulk-model
                        out2bulk <- ode(t = times, y = stateBulk,
                                        func = ModelBulkNutr,
                                        parms = parmsBulk,
                                        rootfun = rootfunBulk,
                                        events = list(func = eventfunBulk,
                                                      root = TRUE,
                                                      terminalroot = 1))
                        if(verbose == TRUE) {
                          print(diagnostics(out2bulk))
                          print(attributes(out2bulk))
                        }
                        EqAfterInvDonorBulk <- tail(out2bulk, 1)[, ]
                        
                        if(plotoutput == 1) {
                          subtitlebulk <- paste0("bR=", bRValue,
                                                 " NI=", NIValue,
                                                 " log10kp=", log10kpValue,
                                                 " log10kn=", log10knValue,
                                                 " e=", NutrConv,
                                                 " w=", wValue,
                                                 " cd=", cdValue,
                                                 " ct=", ctValue,
                                                 " log10gd=", log10gdValue,
                                                 " log10gt=", log10gtValue,
                                                 " gdbulk=", signif(gdbulk, digits = 4),
                                                 " gtbulk=", signif(gtbulk, digits = 4))
                          matplot.deSolve(out2bulk, main = "Bulk-conjugation model", sub = subtitlebulk, ylim = myylim,
                                          log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                                          legend = list(x = "bottomright"))
                          if(verbose == TRUE) {
                            abline(h = extinctionthreshold)
                            abline(v =  attributes(out2bulk)$troot)
                          }
                          grid()
                        }
                      }
                    }
                    
                    MyData[CurrentIteration, c(1:14)] <- unname(c(parmsBulk, Eq))
                    if(runsimulation == 1) {
                      MyData[CurrentIteration, c(15:28)] <- unname(c(
                        EqAfterInvDonor, EqAfterInvDonorBulk))                      
                    } else {
                      MyData[CurrentIteration, c(15:28)] <- rep(NA, 14)
                    }
                    MyData[CurrentIteration, c(29:51)] <- unname(c(
                      EigValEq, DomEigVal, SignDomEigVal, SignEigValEqual,
                      ComplexEigVal, EigValEqBulk, DomEigValBulk,
                      SignDomEigValBulk, SignEigValEqualBulk,
                      ComplexEigValBulk, smallchange, Mytmax, Mytstep))
                  }
                }


BackupMyData <- MyData
colnames(MyData) <- c(
  names(parmsBulk), paste0(names(Eq), "RecipEq"), names(EqAfterInvDonor),
  paste0(names(EqAfterInvDonorBulk), "Bulk"), paste0("EigVal", 1:length(EigValEq)),
  "DomEigVal", "SignDomEigVal", "SignEigValEqual", "ComplexEigVal",
  paste0("EigValBulk", 1:length(EigValEqBulk)),
  "DomEigValBulk", "SignDomEigValBulk", "SignEigValEqualBulk", "ComplexEigValBulk",
  "smallchange", "tmax", "tstep", "PlasmidsInit", "BioInit", "DonorsEq", "PlasmidsEq", "BioEq", "Invasion",
  "DonorsEqBulk", "PlasmidsEqBulk", "BioEqBulk", "InvasionBulk")
if(runsimulation == 1) {
  MyData[, "PlasmidsInit"] <- rep(Dinit, TotalIterations)
  MyData[, "BioInit"] <- MyData[, "RRecipEq"] + MyData[, "PlasmidsInit"]
  MyData[, "DonorsEq"] <- MyData[, "D"] + MyData[, "Mdr"] + MyData[, "Mdt"]
  MyData[, "PlasmidsEq"] <- MyData[, "D"] + MyData[, "Trans"] + MyData[, "Mdr"] +
    2*MyData[, "Mdt"] + MyData[, "Mrt"] + 2*MyData[, "Mtt"]
  MyData[, "BioEq"] <- MyData[, "D"] + MyData[, "R"] + MyData[, "Trans"] + 2*MyData[, "Mdr"] +
    2*MyData[, "Mdt"] + 2*MyData[, "Mrt"] + 2*MyData[, "Mtt"]
  MyData[which(MyData[, "PlasmidsEq"] / MyData[, "BioEq"] > 
                 MyData[, "PlasmidsInit"] / MyData[, "BioInit"]), "Invasion"] <- 1
  MyData[which(MyData[, "PlasmidsEq"] / MyData[, "BioEq"] <= 
                 MyData[, "PlasmidsInit"] / MyData[, "BioInit"]), "Invasion"] <- 0
  MyData[, "DonorsEqBulk"] <- MyData[, "DBulk"]
  MyData[, "PlasmidsEqBulk"] <- MyData[, "DBulk"] + MyData[, "TransBulk"]
  MyData[, "BioEqBulk"] <- MyData[, "DBulk"] + MyData[, "RBulk"] + MyData[, "TransBulk"]
  MyData[which(MyData[, "PlasmidsEqBulk"] / MyData[, "BioEqBulk"] > 
                 MyData[, "PlasmidsInit"] / MyData[, "BioInit"]), "InvasionBulk"] <- 1
  MyData[which(MyData[, "PlasmidsEqBulk"] / MyData[, "BioEqBulk"] <= 
                 MyData[, "PlasmidsInit"] / MyData[, "BioInit"]), "InvasionBulk"] <- 0
} else {
  # Only the parametervalues, plasmid-free equilibrium, and the eigenvalues are stored
  MyData[, "BioInit"] <- MyData[, "RRecipEq"]
  MyData <- MyData[, c(1:14, 29:48, 53)]
}

MyData <- as.data.frame(MyData)
write.csv(MyData, file = paste0(DateTimeStamp, "outputdeSolveChangedContact", ".csv"),
          quote = FALSE, row.names = FALSE)

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = SignDomEigVal)) + 
  ggtitle("Sign dominant eigenvalue of the pair-formation model") +
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0) +
  facet_grid(cd + log10gd ~ ct + log10gt, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputSignDomEigValPair.png"))

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = SignDomEigValBulk)) + 
  ggtitle("Sign dominant eigenvalue of the bulk-conjugation model") +
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0) +
  facet_grid(cd + log10gd ~ ct + log10gt, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputSignDomEigValBulk.png"))

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = SignDomEigVal / SignDomEigValBulk)) + 
  ggtitle("Differences in sign of the dominant eigenvalue") +
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0) +
  facet_grid(cd + log10gd ~ ct + log10gt, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputSignDomEigVals.png"))

summary(MyData$SignDomEigVal - MyData$SignDomEigValBulk)
summary(MyData$SignDomEigVal / MyData$SignDomEigValBulk)

summary(MyData$DomEigVal)
summary(MyData$DomEigValBulk)
DomEigVals <- c(MyData$DomEigVal, MyData$DomEigValBulk)
summary(DomEigVals)

limitseigenvalues <- log10(range(DomEigVals[DomEigVals > 0]))
limitseigenvalues <- c(floor(limitseigenvalues[1]), ceiling(limitseigenvalues[2]))

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(DomEigVal))) + 
  ggtitle("Dominant eigenvalues of the pair-formation model") +
  geom_raster() + 
  scale_fill_gradientn(colours = MyColorBrew, limits = limitseigenvalues) +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputDomEigVal", ".png"))

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(DomEigValBulk))) + 
  ggtitle("Dominant eigenvalues of the bulk-conjugation model") +
  geom_raster() + 
  scale_fill_gradientn(colours = MyColorBrew, limits = limitseigenvalues) +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputDomEigValBulk", ".png"))

limitsbulkrates <- c(floor(min(log10(c(MyData$gtbulk, MyData$gdbulk)))),
                     ceiling(max(log10(c(MyData$gtbulk, MyData$gdbulk)))))

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(gdbulk))) + 
  ggtitle("Bulk-conjugation rate from the donor") +
  geom_raster() + 
  scale_fill_gradientn(colours = MyColorBrew, limits = limitsbulkrates) +
  facet_grid(cd + log10gd ~ ct + log10gt, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputLog10gdbulk.png"))

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(gtbulk))) + 
  ggtitle("Bulk-conjugation rate from the transconjugant") +
  geom_raster() + 
  scale_fill_gradientn(colours = MyColorBrew, limits = limitsbulkrates) +
  facet_grid(cd + log10gd ~ ct + log10gt, labeller = label_both) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)",
       caption = DateTimeStamp)
ggsave(paste0(DateTimeStamp, "outputLog10gtbulk.png"))

## Some controls
any(MyData$time == MyData$tmax) # simulation not complete
length(which(MyData$time == MyData$tmax))
any(MyData$timeBulk == MyData$tmax)
length(which(MyData$timeBulk == MyData$tmax))

any(MyData$R == MyData$BioEq & MyData$SignDomEigVal == 1) # Unstable equilibrium, but invasion not possible
any(MyData$RBulk == MyData$BioEqBulk & MyData$SignDomEigValBulk == 1)
# If the biomass is the same accross simulations, one could use the number of cells instead of fractions of total biomass
summary(MyData$BioEq)
max(MyData$BioEq) / min(MyData$BioEq) 
summary(MyData$BioEqBulk)
max(MyData$BioEqBulk) / min(MyData$BioEqBulk)
summary(MyData$BioEq / MyData$BioEqBulk) # If this differs, comparisons between the 2 models should use fraction, not cell counts
range(MyData$BioEq / MyData$BioEqBulk)

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = BioEq / BioEqBulk)) + 
  ggtitle("Difference in biomass") +
  geom_tile(colour = "white") + 
  scale_fill_gradient2(na.value = "purple") +
  geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
  geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom")

# A <- ggplot(data = MyData, aes(x = log10kp, y = log10(PlasmidsEq))) +
#   ggtitle("Invasion of donor in (R*, Nutr*)") +
#   scale_x_continuous() +
#   scale_y_continuous() +
#   facet_grid(smallchange ~ log10kn, labeller = label_both) +
#   theme(legend.position="bottom") +
#   geom_point(data = MyData, shape = 16, size = 2, aes(color = as.factor(tmax)))
# print(A)
# ggsave(paste0(DateTimeStamp, "outputplotA", ".png"))

###
## Idee: in plaats van runnen met cd = c(0.01, 0.025, 0.05) en ct = c(0.01, 0.025, 0.05)
## en vervolgens met plotten cd ~ ct gebruiken voor facetten, kan ook runnen met
## cd = 0.025, ct = c(0.01, 0.025, 0.05) (dus kleiner dan, gelijk aan, groter dan)
## en zelfdde voor gd en gt, dus gd = 15, gt = 10, 15, 20.

### Plots controleren: facet_grid verschilt
# Note that geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) 
# only takes pair-formation model into account !

### NOTE: ! HARDCODED limits c(-0.3, 0) !
if(length(log10gdSet)*length(log10gdSet) > 1) {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(PlasmidsEq/BioEq))) + 
    ggtitle("Fraction plasmid-bearing bacteria, pair-formation model") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
    geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
    scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(rows = vars(cd, log10gd), cols = vars(ct, log10gt), labeller = label_both) +
    theme(legend.position = "bottom")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionLog.png"))
} else {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(PlasmidsEq/BioEq))) + 
    ggtitle("Fraction plasmid-bearing bacteria, pair-formation model") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50", limits = c(-0.3, 0)) +
    theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
    # geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
    # scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(cd ~ ct, labeller = label_both) +
    theme(legend.position = "bottom") +
    labs(x = "log10(Attachment rate)",
         y = "log10(Detachment rate)")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionLog.png"))
}

ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = PlasmidsEq/BioEq)) + 
  ggtitle("Fraction plasmid-bearing bacteria, pair-formation model") +
  geom_tile(colour = "white") + 
  scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
  theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  # geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
  # scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom") +
  labs(x = "log10(Attachment rate)",
       y = "log10(Detachment rate)")
ggsave(paste0(DateTimeStamp, "outputheatmapFraction.png"))

### NOTE: ! HARDCODED limits c(-0.3, 0) !
if(length(log10gdSet)*length(log10gdSet) > 1) {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(PlasmidsEqBulk/BioEqBulk))) + 
    ggtitle("Fraction plasmid-bearing bacteria, bulk-conjugation model") +
    geom_tile(colour = "white") +
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
    geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
    scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(rows = vars(cd, log10gd), cols = vars(ct, log10gt), labeller = label_both) +
    theme(legend.position = "bottom")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionBulk.png"))
} else {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(PlasmidsEqBulk/BioEqBulk))) + 
    ggtitle("Fraction plasmid-bearing bacteria, bulk-conjugation model") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50", limits = c(-0.3, 0)) +
    theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
    # geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
    # scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(cd ~ ct, labeller = label_both) +
    theme(legend.position = "bottom") +
    labs(x = "log10(Attachment rate)",
         y = "log10(Detachment rate)")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionBulk.png"))
}

min((MyData$PlasmidsEq / MyData$BioEq)[which(MyData$PlasmidsEq / MyData$BioEq != 0)])
max((MyData$PlasmidsEq / MyData$BioEq)[which(MyData$PlasmidsEq / MyData$BioEq != 0)])

min((MyData$PlasmidsEqBulk / MyData$BioEqBulk)[which(MyData$PlasmidsEqBulk / MyData$BioEqBulk != 0)])
max((MyData$PlasmidsEqBulk / MyData$BioEqBulk)[which(MyData$PlasmidsEqBulk / MyData$BioEqBulk != 0)])

sort(unique(((MyData$PlasmidsEq / MyData$BioEq) / (MyData$PlasmidsEqBulk / MyData$BioEqBulk))))

# If the next plot gives problems because of log10(0), could also use log10(... - ...)
# or log10(... + 1)

if(length(log10gdSet)*length(log10gdSet) > 1) {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                            fill = log10((PlasmidsEqBulk/BioEqBulk) / (PlasmidsEq/BioEq)))) + 
    ggtitle("Ratio of the fraction plasmid-bearing bacteria, between the two models") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
    geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
    geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
    scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(rows = vars(cd, log10gd), cols = vars(ct, log10gt), labeller = label_both) +
    theme(legend.position = "bottom")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionPairBulk", ".png"))
} else {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                            fill = log10((PlasmidsEqBulk/BioEqBulk) / (PlasmidsEq/BioEq)))) +
    ggtitle("Ratio of the fraction plasmid-bearing bacteria, between the two models") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
    geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
    geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
    scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(cd ~ ct, labeller = label_both) +
    theme(legend.position = "bottom")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionPairBulk", ".png"))
}

if(length(log10gdSet)*length(log10gdSet) > 1) {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                            fill = log10(PlasmidsEqBulk / PlasmidsEq))) + 
    ggtitle("Ratio of the number of plasmid-bearing bacteria between the two models") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
    geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
    geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
    scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(rows = vars(cd, log10gd), cols = vars(ct, log10gt), labeller = label_both) +
    theme(legend.position = "bottom")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionPairBulk", ".png"))
} else {
  ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                            fill = log10(PlasmidsEqBulk / PlasmidsEq))) +
    ggtitle("Ratio of the number of plasmid-bearing bacteria between the two models") +
    geom_tile(colour = "white") + 
    scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
    geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
    geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
    scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
    facet_grid(cd ~ ct, labeller = label_both) +
    theme(legend.position = "bottom")
  ggsave(paste0(DateTimeStamp, "outputheatmapFractionPairBulk", ".png"))
}

ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                          fill = DonorsEq/BioEq)) +
  ggtitle("Fraction donors at equilibrium, pairmodel") +
  geom_tile(colour = "white") +
  scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
  theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom")
ggsave(paste0(DateTimeStamp, "outputheatmapFractionDonorsPair", ".png"))

ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                          fill = (DonorsEqBulk/BioEqBulk) / (DonorsEq/BioEq))) +
  ggtitle("Ratio of the fraction donors between the two models") +
  geom_tile(colour = "white") +
  scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
  theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom")
ggsave(paste0(DateTimeStamp, "outputheatmapFractionDonorsPairBulk", ".png"))

ggplot(data = MyData, aes(x = log10kp, y = log10kn,
                          fill = (RBulk/BioEqBulk) / (D/BioEq))) +
  ggtitle("Ratio of the fraction donors between the two models") +
  geom_tile(colour = "white") +
  scale_fill_gradientn(colours = MyColorBrew, na.value = "gray50") +
  theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  geom_point(data = MyData[MyData$timeBulk==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom")
ggsave(paste0(DateTimeStamp, "outputheatmapFractionDonorsPairBulk", ".png"))


# Can I set scale limits to compare between plots?

CleanData <- MyData[MyData$time < MyData$tmax, ]
CleanData <- CleanData[CleanData$timeBulk < CleanData$tmax, ]
MyDataBackup <- MyData
MyData <- CleanData

E <- ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = SignEigValEqual)) + 
  ggtitle("pair-formation model") +
  geom_tile(colour = "white") + 
  scale_fill_gradient(low = "red", high = "green") +
  geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ ct, labeller = label_both) +
  theme(legend.position = "bottom")
print(E)
ggsave(paste0(DateTimeStamp, "outputEqualEigenvalues", ".png"), plot = E)

## An alternative way to highlight where equilibrium was not reached, but I could
# not figure out how to include info on the rectangles in the legend, and the
# creation of the raster dataframe requires re-entering the variable names, so
# is more error-prone. On the frames to indicate specific tiles, see
# https://stackoverflow.com/questions/13258454/marking-specific-tiles-in-geom-tile-geom-raster

frames = MyData[MyData$time==MyData$tmax, c("log10kp", "log10kn", "cd")]
ggplot(data = MyData, aes(x = log10kp, y = log10kn, fill = log10(PlasmidsEq))) + 
  geom_rect(colour = "white") + 
  scale_fill_gradient(low = "red", high = "green") +
  geom_rect(data=frames, size=1, fill=NA, colour="black",
            aes(xmin=log10kp - 0.05, xmax=log10kp + 0.05, ymin=log10kn - 0.05, ymax=log10kn + 0.05)) +
  # scale_colour_manual(name = "Ooh look", values = c("black"), labels = "Something cool") +
  facet_grid(cd ~ .) +
  labs(x = "log10kp",
       y = "log10kn",
       title = "pair-formation model", 
       subtitle = "facet by cd ~ .") +
  theme(legend.position = "bottom")

## Andere mogelijk aanpak: subset data in 4 groepen:
## (1) SignDomEig < 0, time == 0 (geen invasie mogelijk, egale kleur want N/A)
## (2) SignDomEig > 0, time >= tmax, invasie mogelijk maar evenwicht niet bereikt,
## dus niet inkleuren maar omkaderen of anderszins markeren
## (3) SignDomEig > 0, time < tmax, PlasmidsEq < (tolerantie?) BioEq (evenwicht
## bereikt, populatie niet volledig plasmiddragend,
## kleurverloop op dit deel van populatie toepassen
## (4) SignDomEig > 0, time < tmax, PlasmidsEq == (tolerantie?) BioEq (evenwicht
## bereikt, populatie volledig plasmiddragend,
## egaal kleuren)


# (1) gaat vanzelf goed via default, want dan geldt PlasmidsEq = 0, dus
# log10(PlasmidsEq) = NA
# (2) is nu geïmplementeerd om de grid cells with a dot te markeren, MAAR DAT
# NOG UIT kleurverloop dataset halen
# (3) is de interessante, ALLEEN DEZE IN KLEURSCHAAL OPNEMEN
# (4) is niet zo interessant, EGAAL KLEUREN EN UIT KLEURSCHAAL HALEN
Pop1 <- MyData[MyData$time == 0, ]
dim(Pop1)
PopRemain <- MyData[MyData$time != 0, ]
dim(PopRemain)

Pop2 <- PopRemain[PopRemain$time >= PopRemain$tmax, ]
dim(Pop2)
PopRemain <- PopRemain[PopRemain$time < PopRemain$tmax, ]
dim(PopRemain)

Pop3 <- PopRemain[PopRemain$PlasmidsEq <= 0.999*PopRemain$BioEq, ]
dim(Pop3)
PopRemain <- PopRemain[PopRemain$PlasmidsEq > 0.999*PopRemain$BioEq, ]
dim(PopRemain)

Pop4 <- PopRemain
dim(Pop4)

## Pop 1 is covered as N/A, but not present in legend
# Pop 2 is covered by geom_point, but should be removed from color
ggplot(data = NULL, aes(x = log10kp, y = log10kn)) + 
  geom_tile(data = MyData, colour = "white") + 
  scale_fill_gradient(low = "red", high = "green") +
  geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ .) +
  labs(x = "log10kp",
       y = "log10kn",
       title = "pair-formation model", 
       subtitle = "facet by cd ~ .") +
  theme(legend.position = "bottom")

# element_text(angle = 45, hjust = 1, vjust = 1, size = 14)
# theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1, vjust = 1))
