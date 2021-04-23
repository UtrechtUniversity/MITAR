################################################################################
## Influence of richness, evenness, diversity, and relatedness on competitive ##
## exclusion of plasmid-free and plasmid-bearing bacteria.                    ##
################################################################################


#### Introduction ####
# How does varying the richness, evenness, and diversity influence the ability
# of bacteria to invade an existing microbiome. How is this altered if the
# invading bacterium carries a plasmid, and conjugation is, or is not,
# relatedness-dependent. We use simulations with a generalised Lotka-Volterra
# model to answer these questions.


#### References ####
# Edelstein-Keshet L. 2005. Mathematical models in biology. Society for
# industrial and applied mathematics.

# Lischke H, Löffler TJ. 2017. Finding all multiple stable fixpoints of n-species
# Lotka-Volterra competition models Theoretical Population Biology 115:24-34.

# May RM. 2001.  Stability and complexity in model ecosystems. Princeton/Oxford:
# Princeton University Press.

# Tokeshi M. 1990. Niche apportionment or random assortment: species abundance
# patterns revisited. Journal of animal ecology 59(3):1129-1146.


#### To do ####

## geteqinfo() ##
# Currently the only the sign and complex part of the largest eigenvalue is used
# to determine the type of equilibrium. However, sometimes the largest
# eigenvalue does not have a complex part when (some of) the other eigenvalues
# do have a complex part. If this affects the equilibrium this should be taken
# into account.
# I check for repeated eigenvalues, but note that the pair a +/- bi are not
# repeated eigenvalues. See p. 133 of Edelstein-Keshet 2005 and section 10.4.3.3 from
# https://eng.libretexts.org/Bookshelves/Industrial_and_Systems_Engineering/Book%3A_Chemical_Process_Dynamics_and_Controls_(Woolf)

## perturbequilibrium() ##
# Abundances frequently grow to infinity because the system does not have an
# inherent carrying capacity. Now I prevent fatal errors during the integration
# by terminating the simulation and indicating infinite growth occurred. Could
# choose to let the object for the output remove before computation, leading to
# NA as final abundance, which will show up in the plot accordingly.


#### Optionally to do ####

## General ##
# I could move the 'niter' argument to the top functions, such that for 1000
# iterations, rnorm is called only once to generate 1000*nspecies growthrates,
# instead of being called 1000 times to generate nspecies growthrates.
# Or use replicate(...) from the apply-family.

# To check if enough iterations are used: run several times for niter iterations,
# if the variation in fraction of stable equilibria is too large, use more
# iterations.

# Create file to store default settings (e.g., sparsity = 0, intdistr =
# selfintdistr = "normal", intsd = selfintsd = 0.1) and write that to a .txt
# file with DateTimeStamp matching to other files, e.g.,
# write.table(x, paste0(DateTimeStamp, "settings"), append = FALSE,
# sep = ",", dec = ".", row.names = TRUE, col.names = TRUE).

# Now I use nested loops, instead I could first create a tibble using
# tidyr::expand_grid(), and then use (l)/(m)apply / purrr:(p)map to iterate over
# all rows?

# Using species as an integer might abolish the need to use as.factor(species)
# when ploting.

## Checking function arguments ##
# See also the remarks on checking function arguments in the 'Optionally to do'
# section below.

# stopifnot() checks conditions of function arguments, but does NOT STOP
# execution of the script if a condition is not met, but issues a warning.
# I should use something from ?conditions that actually terminates execution, or
# use if(!is.numeric(nspecies) || length(nspecies) != 1 || nspecies < 2) {
# stop("'nspecies' must be a length-1 numeric vector with value larger than 1")}
# instead of stopifnot(is.numeric(nspecies), length(nspecies) == 1, nspecies > 1).

# Checking function arguments through stopifnot(...) for every iteration in the
# loop is not needed, and slows down computations (~10%). So instead only use it
# on the function where everything is put together (i.e., the function that is
# eventually called in the main script)? Or create a separate checkarg function
# that has a list containing all variables and criteria, checks the criteria for
# the supplied arguments. Then call that function only once outside the loop? Or
# only check input on the first iteration, by including a checkinput argument
# and use if(checkinput == TRUE) {code to check input}.


## Abundance models ##
# brokenstick()
# diff() works vectorised on matrix columns, but then the breakpoints should be
# sorted in per column before applying diff(), so then still a loop/apply has to
# be used.

# dompreempt()
# Instead of assigning all remaining niche to the last species to obtain the
# user-defined total abundance, Tokeshi 1990 suggests dividing by
# 1 - (1 - 0.75)^nspecies for scaling the geometric series.


## getintmat() ##
# I could use the same interaction matrix for both abundance models.
# See May 2005 on other options for intraspecies interactions (all < 0, all -1, ...).
# Search literature for abundance models applied to a microbiome.

# Add logistic interaction in addition the the linear interaction currently
# implemented (see https://github.com/EgilFischer/FlockMicrobiome for code).


## checkequilibrium() ##
# Instead of a user-defined tmax if showplot = TRUE, I could use a rootfunction
# to stop simulation when the new equilibrium is reached. See ?deSolve::roots


## geteqinfo() ##
# Determining the sign of the real and complex parts can be done vectorised
# outside the loop

## perturbequilibrium() ##
# I could try if supplying the analytic Jacobian to the solver speeds up the
# integration. See Box 1 in Lischke 2017.
   
# Add function to perturb equilibrium when abundances, intmat, ect. are not
# supplied, but instead the parameters to obtain them. Then split this in
# calcStuff function (to be used also in main script [with (l/m)apply?]) and
# PerturbFunction ?



#### Loading required libraries ####
library(deSolve)   # checkequilibrium calls ode() if showplot == TRUE
library(dplyr)     # checkequilibrium and perturbequilibrium call near()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()


#### Settings and defining parameterspace ####

# Simulation settings
niter <- 100
saveplots <- TRUE
smallstate <- 1e-20 # States are set to 0 if they become smaller than smallstate
smallchange <- 1e-10 # If the sum of absolute rates of change is equal to
# smallchange, equilibrium is assumed to be reached and integration is terminated

# Define parameter space
totalabun <- 1
nspeciesset <- c(2, 4, 6)
abunmodelset <- c("brokenstick", "dompreempt")
intmeanset <- seq(from = -0.8, to = 0.8, by = 0.1)
selfintmeanset <- seq(from = -0.8, to = -0.5, by = 0.025)
costset <- c(0.01, 0.20)
conjrateset <- c(0.01, 0.05, 0.1)
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")

# Settings for testing code
niter <- 2
saveplots <- TRUE
totalabun <- 1
nspeciesset <- c(2, 4)
abunmodelset <- c("brokenstick", "dompreempt")
intmeanset <- seq(from = -0.8, to = 0.8, by = 0.8)
selfintmeanset <- seq(from = -0.8, to = 0.8, by = 0.8)
costset <- c(0.01, 0.20)
conjrateset <- c(0.01, 0.1)
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")


#### Functions ####

# The generalised Lotka-Volterra model. n is a vector of species densities (cell
# mL^-1), dn/dt is a vector of the change in these densities per time (cell
# mL^-1 time^-1), growthrate is a vector of growth rates (h^-1), intmat is a
# matrix of scaled interaction coefficients with units mL cell^-1 h^-1, where
# element aij represents cij * ri / Ki in the textbook-notation given below,
# such that elements aii on the diagonal equal ri / Ki, with Ki being the
# carrying capacity of species i.

# Textbooks (e.g., Eq. 9 in Edelstein-Keshet 2005, p. 224) explicitly include
# the carrying capacities into gLV models:
# dn1/dt = r1*n1*((K1 -     n1 - c12*n2)/K1) = r1*n1*(1 - (    n1 + c12*n2)/K1)
# dn2/dt = r2*n2*((K2 - c21*n1 -     n2)/K2) = r2*n2*(1 - (c21*n1 +     n2)/K2)
# In matrix-notation this becomes: dn/dt <- r*n*(1 - (1/K) * intmat %*% n).
# Ki is the carrying capacity of species 1 (cells mL^-1), and the interaction
# matrix gives dimensionless interaction coefficients where element cij gives
# the decline (if cij is positive) or increase (if cij is negative) in the
# growth rate of species i caused by one individual of species j. The diagonal
# entries of intmat are the intraspecies interaction coefficients cii and should
# be -1 to ensure that species in isolation grow to their carrying capacities K.

# Instead of normal multiplication using n, matrix multiplication using
# diag(n) can be used. This results in diag(n) %*% (growthrates + intmat %*% n).
gLV <- function(t, n, parms) {
  with(parms, {
    dn <- n*(growthrate + intmat %*% n)
    return(list(c(dn)))
  })
}

# Define the generalised Lotka-Volterra model with conjugation. Vector n with
# species abundances is split into vector S0 for plasmid-free bacteria and
# vector S1 for plasmid-bearing bacteria. cost is a vector giving the absolute
# reduction in growth rate when carrying a plasmid. conjmat is a matrix with
# conjugation rates, where element cij gives the rate of conjugation from
# plasmid-bearing species j to plasmid-free species i.
gLVConj <- function(t, n, parms) {
  with(parms, {
  S0 <- n[1:nspecies]
  S1 <- n[(nspecies+1):(2*nspecies)]
  
  dS0 <- S0*(growthrate        + intmat %*% (S0 + S1)) - (conjmat %*% S1) * S0
  dS1 <- S1*(growthrate - cost + intmat %*% (S0 + S1)) + (conjmat %*% S1) * S0
  
  dn <- c(dS0, dS1)
  return(list(dn))
  })
}

# Calculate absolute species abundances from the specified total abundance if
# species abundances are proportional to the length of fragments of a stick that
# is broken randomly at nspecies - 1 points, following the broken stick model
# (= MacArthur faction model) in the description of Tokeshi (1990).
brokenstick <- function(nspecies, totalabun, takelimit = TRUE) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
           length(totalabun) == 1, totalabun > 0)
  if(takelimit == TRUE) {
    niterabun <- 1e3
    abunmat <- matrix(data = NA, nrow = niterabun, ncol = nspecies)
    for(iterindex in 1:niterabun) {
      # sorting to get positive differences in the next step
      breakpoints <- sort(runif(nspecies - 1, min = 0, max = totalabun))
      # sorting because otherwise taking the mean does not make sense
      abunmat[iterindex, ] <- sort(diff(c(0, breakpoints, totalabun)))
    }
    abun <- sort(colMeans(abunmat), decreasing = TRUE)
  } else {
    breakpoints <- sort(runif(nspecies - 1, min = 0, max = totalabun))
    abun <- sort(diff(c(0, breakpoints, totalabun)), decreasing = TRUE)
  }
  return(abun)
}

# Calculate abundances following the dominance preemption model. The first
# species occupies (preempts) more than half of the total niche, and each
# subsequent species occupies more than half of the remainder. The last species
# is assigned all of the niche that remains, to obtain the user-defined total
# abundance. Over many iterations, each species preempts on average (0.5 + 1)/2
# = 0.75 of the remainder, such that the model converges to the geometric series
# with k = 0.75 for all but the last species. This geometric model is used
# instead of the dominance preemption model with the default takelimit = TRUE.
dompreempt <- function(nspecies, totalabun, takelimit = TRUE) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
           length(totalabun) == 1, totalabun > 0)
  abun <- rep(NA, nspecies)
  remainingabun <- totalabun  
  if(takelimit == TRUE) {
    for(speciesindex in 1:nspecies) {
      abun.temp <- 0.75*remainingabun
      abun[speciesindex] <- abun.temp
      remainingabun <- remainingabun - abun.temp
    }
    abun[speciesindex] <- abun[speciesindex] + remainingabun
  } else {
    for(speciesindex in 1:nspecies) {
      abun.temp <- runif(1, min = 0.5, max = 1)*remainingabun
      abun[speciesindex] <- abun.temp
      remainingabun <- remainingabun - abun.temp
    }
    abun[speciesindex] <- abun[speciesindex] + remainingabun
  }
  return(abun)
}

# Create a matrix of scaled interaction coefficients for nspecies species with
# units mL cell^-1 h^-1, where element aij represents cij * ri / Ki in the
# textbook-notation (see comments on the generalised Lotka-Volterra model given
# above), such that elements aii on the diagonal equal ri / Ki, with Ki being
# the carrying capacity of species i. The fraction of sparse interspecies
# interactions can be set through 'sparsity', with sparsity = 0 leading to a
# fully connected matrix, and sparsity = 1 leading to all off-diagonal entries
# equal to 0. Off-diagonal entries are interspecies interaction coefficients
# drawn from the distribution 'intdistr'. Diagonal entries are self-interactions
# drawn from the distribution 'selfintdistr', which is truncated to obtain
# negative self-interactions. To get fixed values for the interactions, choose
# the uniform distribution and provide the desired value both as the minimum and
# maximum of the range. The other arguments specify the distributions from which
# interaction coefficients are drawn.
getintmat <- function(nspecies, sparsity = 0,
                      intdistr = "normal", intmean = 0, intsd = 0.4,
                      intrange = c(-0.8, 0.8),
                      selfintdistr = "normal", selfintmean = -0.65, selfintsd = 0.075,
                      selfintrange = c(-0.8, -0.5)) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
           is.numeric(sparsity), length(sparsity) == 1,
           sparsity >= 0, sparsity <= 1)
  switch(intdistr,
         normal = {
           stopifnot(length(intmean) == 1, length(intsd) == 1, intsd > 0)
           intcoefs <- rnorm(n = nspecies^2, mean = intmean, sd = intsd)
         },
         uniform = {
           stopifnot(length(intrange) == 2, intrange[1] <= intrange[2])
           intcoefs <- runif(n = nspecies^2, min = intrange[1],
                             max = intrange[2])
         },
         {
           warning("'intdistr' should be 'normal' or 'uniform'.")
           intcoefs <- NULL
         }
  )
  intmat <- matrix(intcoefs, nrow = nspecies, ncol = nspecies)
  
  switch(selfintdistr,
         normal = {
           stopifnot(length(selfintmean) == 1, length(selfintsd) == 1,
                     selfintsd > 0)
           diag(intmat) <- rnorm(n = nspecies, mean = selfintmean, sd = selfintsd)
           # Draw variates from a truncated normal distribution to ensure
           # negative self-interactions without having to redraw for positive
           # deviates, using rtnorm() from the package TruncatedNormal.
           rtnorm(n = nspecies, mu = selfintmean, sd = selfintsd,
                  lb = -Inf, ub = 0, method = "invtransfo")
           },
         uniform = {
           stopifnot(length(selfintrange) == 2, selfintrange[1] <= selfintrange[2],
                     selfintrange[2] < 0)
           diag(intmat) <- runif(n = nspecies,
                                 min = selfintrange[1], max = selfintrange[2])
         },
         {
           warning("'selfintdistr' should be 'normal' or 'uniform'.")
           diag(intmat) <- NULL
         }
  )
  
  if(sparsity > 0) {
    nsparseint <- round(sparsity*(nspecies^2 - nspecies))
    # Create indexmatrix with column- and row indices
    indexmat <- matrix(c(rep(1:nspecies, nspecies),
                         rep(1:nspecies, each = nspecies)),
                       ncol = 2, dimnames = list(NULL, c("row", "column")))
    # Only keep rows of indexmat specifying off-diagonal entries of intmat, to
    # ensure that self-interaction coefficients never become sparse.
    indexmat <- indexmat[indexmat[, "row"] != indexmat[, "column"], ]
    # Sample rows of indexmat to get index of matrix entries that become sparse
    sparse.index <- sample(1:(dim(indexmat)[1]), nsparseint)
    # Only keep rows of indexmat that were drawn from the sample
    indexmat <- indexmat[sparse.index, ]
    # Set entries of the interaction matrix indicated by indexmat to zero.
    intmat[indexmat] <- 0
  }
  return(intmat)
}

# Calculate the required growth rates to obtain the specified species abundances
# at equilibrium given the interaction matrix. The growth rate can be negative
# if growth is assumed to be slower than washout.
getgrowthrate <- function(abundance, intmat) {
  stopifnot(length(abundance) == dim(intmat)[1])
  growthrate <- -intmat %*% abundance
  if(any(is.complex(growthrate))) {
    countcomplex <- length(which(Im(growthrate) != 0))
    warntext <- paste(countcomplex, "growth rates contain an imaginary part.")
    warning(warntext)
  }
  return(growthrate)
}

# Check if the analytically identified equilibrium is indeed an equilibrium, by
# checking if the derivative at the presumed equilibrium is zero. If showplot =
# TRUE, a time course starting from the presumed equilibrium is shown.
checkequilibrium <- function(abundance, intmat, growthrate,
                             printderivatives = FALSE,
                             showplot = FALSE, tmax = 100, tstep = 0.1) {
  derivatives <- unlist(gLV(t = 0, n = abundance,
                            parms = list(growthrate = growthrate, intmat = intmat)))
  atequilibrium <- all(near(0, derivatives))
  if(printderivatives == TRUE) {
    print(paste("Derivatives:", paste0(signif(derivatives, 4), collapse = ", ")),
          quote = FALSE)
  }
  if(showplot == TRUE) {
    times <- seq(from = 0, to = tmax, by = tstep)
    out <- ode(y = abundance, t = times, func = gLV,
               parms = list(growthrate = growthrate, intmat = intmat))
    ylim <- c(0, 1.1*max(out[, -1]))
    matplot.deSolve(out, ylim = ylim, lwd = 2,
                    lty = 1, ylab = "Abundance")
    grid()
  }
  return(atequilibrium)
}

getconjmat <- function(nspecies, conjrate) {
  conjmat <- matrix(rep(conjrate, nspecies^2),
                    nrow = nspecies, ncol = nspecies)
}

# Get largest eigenvalues of the model without and with conjugation
geteqinfo <- function(abundance, intmat,
                      growthrate, cost, conjmat) {
  
  eigval <- eigen(x = jacobian.full(y = abundance, func = gLV,
                      parms = list(growthrate = growthrate, intmat = intmat)),
                  symmetric = FALSE, only.values = TRUE)$values
  eigvalRep <- any(duplicated(eigval))
  # Using sort(eigval) to get eigenvalue with largest real part, because
  # max(eigval) does not work if eigval is complex. Complex values are sorted
  # first by the real part, then the imaginary part.
  eigval <- sort(eigval, decreasing = TRUE)[1]
  eigvalRe <- Re(eigval)
  eigvalIm <- Im(eigval)
  eigvalReSign <- sign(eigvalRe)
  eigvalImSign <- sign(eigvalIm)
  
  eigvalconj <- eigen(x = jacobian.full(y = c(abundance, rep(0, nspecies)),
                                        func = gLVConj,
                                        parms = list(growthrate = growthrate, intmat = intmat,
                                                     cost = cost, conjmat = conjmat)),
                      symmetric = FALSE, only.values = TRUE)$values
  eigvalconjRep <- any(duplicated(eigvalconj))
  eigvalconj <- sort(eigvalconj, decreasing = TRUE)[1]
  eigvalconjRe <- Re(eigvalconj)
  eigvalconjIm <- Im(eigvalconj)
  eigvalconjReSign <- sign(eigvalconjRe)
  eigvalconjImSign <- sign(eigvalconjIm)
  
  eqinfo <- c(eigvalRe, eigvalIm, eigvalReSign, eigvalImSign, eigvalRep,
              eigvalconjRe, eigvalconjIm, eigvalconjReSign, eigvalconjImSign, eigvalconjRep)
  return(eqinfo)
}

# Root-functions and event-functions
# A value returned by the rootfunction becomes zero ('a root is found') in three
# situations: (1) if the sum of absolute rates of change is equal to threshold
# smallchange (indicating that equilibrium has been reached), (2) if abundances
# get very large (indicating unrestricted growth, which would eventually lead to
# errors during integration), (3) if any of the state variables gets smaller
# than smallstate.
# When a root is found, the event-function is called to set any state variables
# smaller than smallstate to 0.
# In the first two cases, the simulation should terminate (in the first case
# because equilibrium has been reached, in the second case to prevent errors
# during integration because states get to +Inf or stepsize gets to 0). This is
# achieved by specifying terminalroot = c(1, 2) within ode(...). See
# help(events) and help(lsodar) (both in the deSolve package) for background
# information and examples.
# COULD replace state[state <= smallstate] in eventfun with state[state <= smallstate]
# to prevent constantly re-copying the zeros ?
rootfun <- function(t, state, p) {
  c(sum(abs(unlist(gLV(t, state, p)))) - smallchange,
    sum(state) - 1e10*totalabun,
    state - smallstate)
}

rootfunconj <- function(t, state, p) {
  c(sum(abs(unlist(gLVConj(t, state, p)))) - smallchange,
    sum(state) - 1e10*totalabun,
    state - smallstate)
}

eventfun <- function(t, state, p) {
  state[state <= smallstate] <- 0
  return(state)
}

# Perturb equilibrium by adding a small number of bacteria to an equilibrium.
# Input:
# abundance: a numeric vector of species abundances at equilibrium (a warning is
#   issued if they are not at equilibrium). If the model 'gLVConj' is selected,
#   the abundances of the plasmid-bearing populations should be included in the
#   vector. If they are not 0, a warning is issued.
# intmat, growthrate, and conjmat are supposed to be the output of getintmat(),
#   getgrowthrate(), and getconjmat()
# cost is a vector of plasmid costs in growth rate. 
# model should be 'gLV' (no plasmids modelled) or 'gLVConj' (to include plasmid-
#   bearing populations in the model).
# pertpop is a character vector with the name(s) of the population(s) to be
#   perturbed. Acceptable names are those of specific populations such as "R1"
#   or "P2", and the groups of populations 'all', 'R' and 'P' to denote all, all
#   plasmid-free, and all plasmid-bearing populations, respectively. These can
#   be combined, for example using pertpop = c("R", "P1")
# pertmagn gives the absolute increase in populations for the pertubation
# tmax and tstep give the timesteps at which abundances should be calculated
#   (since variable step-size methods are used, those are not the only times
#   that integration occurs)
# If showplot == TRUE, the result is plotted, which is slow
# If verbose is TRUE, abundances before and after perturbation, and their
#   differences, are printed.
# Abundances grow to infinity when positive feedback is present, because the
#   system does not have an inherent carrying capacity. To prevent fatal errors
#   during the integration caused by state variables going to +Inf, or stepsizes
#   to 0, a rootfunction is used to terminate the integration if abundances get
#   very large. A warning is issued and the variable 'infgrowth' is set to 1 if
#   this occurs.
# 
# Returns abunfinal, a list containing:
#   - abunfinalR with abundances for plasmid-free species,
#   - abunfinalP with abundances for plasmid-bearing species (if model ==
#     "gLVConj"), or NULL (if model == "gLV")
#   - infgrowth indicating if infinite growth was detected
# result1 <- perturbequilibrium(abundance, intmat, growthrate, cost, conjmat, "gLV", "R1")
perturbequilibrium <- function(abundance, intmat, growthrate, cost, conjmat,
                             model, pertpop, pertmagn = 1e-6,
                             tmax = 100, tstep = 0.1, showplot = TRUE,
                             verbose = TRUE, suppresswarninfgrowth = FALSE) {
  
  # Name abundances, set line type and colors, get derivatives of initial state
  # in the plasmid-free model.
  if(model == "gLV") {
    nspecies <- length(abundance)
    names(abundance) <- paste0(rep("R", nspecies), 1:nspecies)
    lty <- 1
    col <- mycol[1:nspecies]
    derivatives <- unlist(
      gLV(n = abundance, parms = list(growthrate = growthrate, intmat = intmat))
    )
  }
  
  if(model == "gLVConj") {
    nspecies <- length(abundance)/2
    names(abundance) <- c(paste0(rep("R", nspecies), 1:nspecies),
                          paste0(rep("P", nspecies), 1:nspecies))
    lty <- rep(c(1, 2), each = nspecies)
    col <- rep(mycol[1:nspecies], 2)
    if(!all(near(abundance[(nspecies + 1):(2*nspecies)], 0))) {
      message("Initial state is NOT plasmid-free.")
    }
    derivatives <- unlist(
      gLV(n = abundance[1:nspecies],
          parms = list(growthrate = growthrate, intmat = intmat)
      )
    )
  }
  
  # Warn if initial plasmid-free state is not an equilibrium in the model
  # without plasmids
  if(!all(near(derivatives, 0))) {
    warntext <- paste("Initial (unperturbed) state is NOT an equilibrium in the
                      plasmid-free model! Derivatives are:",
                      paste(signif(derivatives), collapse = ", ")
    )
    warning(warntext)
  }
  
  # Get names of populations to perturb
  if("all" %in% pertpop) {
    pertpop <- names(abundance)
  }
  if("R" %in% pertpop) {
    pertpop <- pertpop[-which(pertpop == "R")]
    pertpop <- unique(c(pertpop, paste0(rep("R", nspecies), 1:nspecies)))
  }
  if("P" %in% pertpop) {
    pertpop <- pertpop[-which(pertpop == "P")]
    pertpop <- unique(c(pertpop, paste0(rep("P", nspecies), 1:nspecies)))
  }
  
  # Only perturb populations that exist
  if(!all(pertpop %in% names(abundance))) { 
    warning("Neglecting non-existent population(s) specified to perturb!")
    pertpop <- pertpop[which(pertpop %in% names(abundance))]
  }
  
  # Create perturbed abundances
  abuninit <- abundance
  abunpert <- abundance
  abunpert[pertpop] <- abunpert[pertpop] + pertmagn
  
  if(verbose == TRUE) {
    print("abunpert =")
    print(abunpert)
  }
  
  # Perturb equilibrium
  times <- seq(from = 0, to = tmax, by = tstep)
  if(model == "gLV") {
    out <- ode(y = abunpert, t = times, func = gLV,
               parms = list(growthrate = growthrate, intmat = intmat),
               rootfun = rootfun,
               events = list(func = eventfun, root = TRUE, terminalroot = c(1, 2)))
  }
  if(model == "gLVConj") {
    out <- ode(y = abunpert, t = times, func = gLVConj,
               parms = list(growthrate = growthrate, intmat = intmat,
                            cost = cost, conjmat = conjmat),
               rootfun = rootfunconj,
               events = list(func = eventfun, root = TRUE, terminalroot = c(1, 2)))
  }
  abunfinal <- tail(out, 1)[, -1]
  names(abunfinal) <- names(abundance)
  
  # Assume infinite growth occurred if a root was triggered because abundances
  # became very large
  infgrowth <- 0
  eqreached <- 0
  if(!is.null(attributes(out)$troot)) {
    # One or more roots found
    
    if(any(attributes(out)$indroot == 1)) {
      eqreached <- 1  
    }
    
    if(any(attributes(out)$indroot == 2)) {
      infgrowth <- 1
      if(suppresswarninfgrowth != TRUE) {
        warning(
          paste0("Integration was terminated at time = ",
                round(attributes(out)$troot[which(attributes(out)$indroot == 2)], 2),
                ", when the sum of abundances became",
                "\nlarger than 1e10 times the initial total abundance, indicating unbounded growth.",
                "\nAbundances then were ",
                paste(names(abunfinal), "=", signif(abunfinal), collapse = ", "))
        )
      }
    }
  }
  
  if(eqreached == 0) {
    warning("Equilibrium has not been reached. Final abundances were\n",
            paste(names(abunfinal), "=", abunfinal, collapse = ", "),
            "Increase tmax to prevent this?")
  }
  
  if(showplot == TRUE) {
    subtitle <- paste0(abunmodel, ", intmean=", intmean, ", selfintmean=", selfintmean, ", cost=",
                       cost, ", conjrate=", conjrate)
    matplot.deSolve(out, lty = lty, col = col, ylab = "Abundance",
                    log = "y", sub = subtitle, lwd = 2, legend = list(x = "bottomright"))
    grid()
    abline(h = abuninit)
  }
  
  if(verbose == TRUE) {
    print(paste("troot=", attributes(out)$troot))
    print("iroot=")
    print(attributes(out)$indroot)
    print(paste("infgrowth=", infgrowth))
    print(paste("eqreached=", eqreached))
    
    abschange <- abunfinal - abuninit
    relchange <- abunfinal / abuninit
    print("Species abundances, and their changes:", quote = FALSE)
    print(rbind(tested = abuninit, at_perturbation = abunpert, final = abunfinal,
                absolute_change = abschange, relative_change = relchange))
  }
  
  if(model == "gLV") {
    abunfinal <- list(abunfinalR = abunfinal,
                      abunfinalP = NULL,
                      infgrowth = infgrowth, eqreached = eqreached)
  } else {
    abunfinal <- list(abunfinalR = abunfinal[1:nspecies],
                      abunfinalP = abunfinal[(nspecies + 1):(2*nspecies)],
                      infgrowth = infgrowth, eqreached = eqreached)
  }
  return(abunfinal)
}

# Function to create plots
CreatePlot <- function(dataplot = plotdata, xvar = "intmean", yvar = "selfintmean",
                       fillvar, filltitle, filltype = "discrete", limits = NULL, 
                       labx = "Mean interaction coefficient",
                       laby = "Mean selfinteraction coefficient",
                       mytag = NULL, addstamp = FALSE, diagional = "none",
                       facetx = "modelcode", facety = "nspecies", as.table = TRUE,
                       marginx = NULL, marginy = NULL, base_size = 11,
                       rotate_legend = FALSE,
                       save = saveplots, filename = NULL, ...) {
  caption <- paste(unique(dataplot$niter), "iterations")
  if(exists("DateTimeStamp") == FALSE) {
    DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
    if(addstamp == TRUE) {
      warning("DateTimeStamp created to include in plot does not correspond to filename of the dataset")
    }
  }
  if(addstamp == TRUE) {
    caption <- paste0(caption, ", ", DateTimeStamp)
  }
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar),
              subtitle = subtitle) + 
    theme_bw(base_size = base_size) +
    geom_raster() +
    scale_x_continuous() +
    scale_y_continuous() +
    coord_fixed(ratio = 1, expand = FALSE) +
    facet_grid(as.formula(paste(facety, "~", facetx)), as.table = as.table,
               labeller = mylabeller) +
    theme(legend.position = "bottom") +
    labs(x = labx, y = laby, caption = caption, tag = mytag)
  if(!is.null(marginx)) {
    p <- p + theme(strip.text.x = element_text(margin = margin(marginx)))
  }
  if(!is.null(marginy)) {
    p <- p + theme(strip.text.y = element_text(margin = margin(marginy)))
  }
  if(filltype == "discrete") {
    p <- p + scale_fill_viridis_d(filltitle, limits = if(is.null(limits)) {
      as.factor(c(-1, 1))
    } else {
      factor(limits)}, labels = filllabels)
  }
  if(filltype == "binned") {
    p <- p + scale_fill_viridis_b(filltitle, breaks = limits)
  }
  if(filltype == "continuous") {
    p <- p + scale_fill_viridis_c(filltitle, limits = limits)
  }
  if(rotate_legend == TRUE) {
    p <- p + guides(fill = guide_colourbar(label.hjust = 0.4, label.vjust = 0.5,
                                           label.theme = element_text(angle = 90)))
  }
  if(diagional == "both" | diagional == "major") {
    p <- p + geom_abline(intercept = 0, slope = -1, col = "white", size = 1.1)
  }
  if(diagional == "both" | diagional == "minor") {
    p <- p + geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1)
  }
  print(p)
  if(save == TRUE) {
    if(is.null(filename)) {
      filename <- gsub("/", ".", fillvar)
      filename <- gsub(" ", "", filename)
      filename <- paste0(filename, filltype)
    }
    filename <- paste0(DateTimeStamp, filename, ".png")
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      ggsave(filename)
    }
  }
}


#### Testing functions ####
# 
# nspecies <- 4
# abunbrokenstick <- brokenstick(nspecies = nspecies, totalabun = totalabun, takelimit = TRUE)
# abundompreempt <- dompreempt(nspecies = nspecies, totalabun = totalabun, takelimit = TRUE)
# abunbrokenstick
# abundompreempt
# 
# (intmat1 <- getintmat(nspecies = nspecies))
# (growthratebrokenstick <- getgrowthrate(abundance = abunbrokenstick, intmat = intmat1))
# (growthratedompreempt <- getgrowthrate(abundance = abundompreempt, intmat = intmat1))
# checkequilibrium(abundance = abunbrokenstick, intmat = intmat1,
#                  growthrate = growthratebrokenstick, printderivatives = TRUE, showplot = TRUE) # At equilibrium
# checkequilibrium(abundance = 0.9*abunbrokenstick, intmat = intmat1,
#                  growthrate = growthratebrokenstick, printderivatives = TRUE, showplot = TRUE) # Not at equilibrium
# 
# # geteqinfo returns eigvalRe, eigvalIm, eigvalReSign, eigvalImSign, and eigvalRep
# geteqinfo(abundance = abunbrokenstick, intmat = intmat1, growthrate = growthratebrokenstick)
# geteqinfo(abundance = abundompreempt, intmat = intmat1, growthrate = growthratedompreempt)


#### Running the simulations ####
set.seed(seed = 314, kind = "default", normal.kind = "default", sample.kind = "default")

# Create matrix to store data
nrowplotdata <- length(nspeciesset)*length(abunmodelset)*
  length(intmeanset)*length(selfintmeanset)*length(costset)*length(conjrateset)
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = 32)
colnames(plotdata) <- c("niter", "nspecies", "modelcode",
                       "intmean", "selfintmean", "cost", "conjrate",
                       "mingrowthrate", "meangrowthrate", "maxgrowthrate",
                       "fracstable", "fracreal", "fracrep",
                       "fracstableconj", "fracrealconj", "fracrepconj",
                       "fracinfgrowth", "fracinfgrowthconj",
                       "fraceqreached", "fraceqreachedconj",
                       "minR", "meanR", "medianR", "maxR",
                       "minRconj", "meanRconj", "medianRconj", "maxRconj",
                       "minPconj", "meanPconj", "medianPconj", "maxPconj")

nrowdatatotal <- length(abunmodelset)*length(intmeanset)*
  length(selfintmeanset)*length(costset)*length(conjrateset)*niter*
  sum(nspeciesset)

mydatatotal <- matrix(data = NA, nrow = nrowdatatotal, ncol = 29)
indexmydatatotal <- 1
maxnspecies <- max(nspeciesset)

# system.time({
# Run simulations
rowindexplotdata <- 1
rowindexmydata <- 1
for(nspecies in nspeciesset) {
  for(conjrate in conjrateset) {
  conjmat <- getconjmat(nspecies = nspecies,
                        conjrate = conjrate)
  
  for(abunmodel in abunmodelset) {

    if(abunmodel == "brokenstick") {
      abundance <- brokenstick(nspecies = nspecies, totalabun = totalabun,
                               takelimit = TRUE)
      modelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt(nspecies = nspecies, totalabun = totalabun,
                              takelimit = TRUE)
      modelcode <- 2
    }
    
    for(intmean in intmeanset) {
      print(paste0("nspecies = ", nspecies, ", conjrate = ", conjrate,
                   ", abundance model = ", abunmodel, ", intmean = ", intmean,
                   ": started at ", Sys.time()), quote = FALSE)
      for(selfintmean in selfintmeanset) {
        for(cost in costset) {
        nrowmydata <- niter * nspecies
        mydata <- matrix(data = NA, nrow = nrowmydata, ncol = 29)
        for(iter in 1:niter) {
          intmat <- getintmat(nspecies = nspecies,
                              intmean = intmean, selfintmean = selfintmean)
          
          growthrate <- getgrowthrate(abundance = abundance, intmat = intmat)

          eqinfo <- geteqinfo(abundance = abundance, intmat = intmat,
                              growthrate = growthrate, cost = cost,
                              conjmat = conjmat)
          
          abunfinal <- perturbequilibrium(abundance = abundance, intmat = intmat,
                                        growthrate = growthrate, cost = cost,
                                        conjmat = conjmat,
                                        model = "gLV", pertpop = "all", tmax = 1e3,
                                        showplot = FALSE, verbose = FALSE)
          infgrowth <- abunfinal$infgrowth
          eqreached <- abunfinal$eqreached
          
          if(infgrowth == 0) {
            abunR <- sum(abunfinal$abunfinalR)
          } else {
            abunR <- NA
          }
          
          # To simulate invasion of the plasmid-free equilibrium with plasmids,
          # the abundances of the plasmid-free populations have to be appended
          # to the abundances of the plasmid-bearing populations
          abunfinalconj <- perturbequilibrium(
            abundance = c(abundance, rep(0, nspecies)),
            intmat = intmat,
            growthrate = growthrate, cost = cost,
            conjmat = conjmat,
            model = "gLVConj", pertpop = "P", tmax = 1e3,
            showplot = FALSE, verbose = FALSE)
          infgrowthconj <- abunfinalconj$infgrowth
          eqreachedconj <- abunfinalconj$eqreached
          
          # It does not make sense to store abundances for infinite growth, so
          # record those as NA
          if(infgrowthconj == 0) {
            abunRconj <- sum(abunfinalconj$abunfinalR)
            abunPconj <- sum(abunfinalconj$abunfinalP)
          } else {
            abunRconj <- NA
            abunPconj <- NA
          }
          mydata[(1 + nspecies*(iter - 1)):(nspecies*iter), ] <- cbind(
            rep(niter, nspecies),
            rep(nspecies, nspecies),
            rep(modelcode, nspecies),
            rep(intmean, nspecies),
            rep(selfintmean, nspecies),
            rep(cost, nspecies),
            rep(conjrate, nspecies),
            rep(iter, nspecies),
            1:nspecies,
            abundance,
            diag(intmat),
            growthrate,
            matrix(rep(eqinfo, nspecies), nrow = nspecies, byrow = TRUE),
            rep(infgrowth, nspecies),
            rep(infgrowthconj, nspecies),
            rep(eqreached, nspecies),
            rep(eqreachedconj, nspecies),
            rep(abunR, nspecies),
            rep(abunRconj, nspecies),
            rep(abunPconj, nspecies)
          )
        }
        mydatatotal[indexmydatatotal:(indexmydatatotal + nrowmydata - 1), ] <- mydata
        indexmydatatotal <- indexmydatatotal + nrowmydata
        
        colnames(mydata) <- c("niter", "nspecies", "modelcode",
                              "intmean", "selfintmean", "cost", "conjrate",
                              "iter", "species", "abundance",
                              "selfint", "growthrate",
                              "eigvalRe", "eigvalIm",
                              "eigvalReSign", "eigvalImSign", "eigvalRep",
                              "eigvalconjRe", "eigvalconjIm",
                              "eigvalconjReSign", "eigvalconjImSign", "eigvalconjRep",
                              "infgrowth", "infgrowthconj",
                              "eqreached", "eqreachedconj",
                              "abunR", "abunRconj", "abunPconj")
        
        # Get proportions of stable and non-oscillating equilibria, and repeated eigenvalues
        fracstable <- mean(mydata[, "eigvalRe"] < 0)
        fracreal <- mean(mydata[, "eigvalImSign"] == 0)
        fracrep <- mean(mydata[, "eigvalRep"] != 0)
        fracstableconj <- mean(mydata[, "eigvalconjRe"] < 0)
        fracrealconj <- mean(mydata[, "eigvalconjImSign"] == 0)
        fracrepconj <- mean(mydata[, "eigvalconjRep"] != 0)
        fracinfgrowth <- mean(mydata[, "infgrowth"] != 0)
        fracinfgrowthconj <- mean(mydata[, "infgrowthconj"] != 0)
        fraceqreached <- mean(mydata[, "eqreached"] != 0)
        fraceqreachedconj <- mean(mydata[, "eqreachedconj"] != 0)
        
        summaryabunR <- summary(mydata[, "abunR"])[c("Min.", "Median", "Mean", "Max.")]
        summaryabunRconj <- summary(mydata[, "abunRconj"])[c("Min.", "Median", "Mean", "Max.")]
        summaryabunPconj <- summary(mydata[, "abunPconj"])[c("Min.", "Median", "Mean", "Max.")]
        
        plotdata[rowindexplotdata, ] <- c(niter, nspecies, modelcode,
                                          intmean, selfintmean, cost, conjrate,
                                          min(mydata[, "growthrate"]),
                                          mean(mydata[, "growthrate"]),
                                          max(mydata[, "growthrate"]),
                                          fracstable, fracreal, fracrep,
                                          fracstableconj, fracrealconj, fracrepconj,
                                          fracinfgrowth, fracinfgrowthconj,
                                          fraceqreached, fraceqreachedconj,
                                          summaryabunR, summaryabunRconj,
                                          summaryabunPconj)
        rowindexplotdata <- rowindexplotdata + 1
        }
      }
    }
  }
  }
}
# })
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)
colnames(mydatatotal) <- colnames(mydata)

#### Reading previously saved data from a .csv-file ####
## To read data from csv-file, uncomment this section and fill in the 
# needed datetimestamp
# filename <- "2021_04_16_17_13multispecies.csv"
# plotdata <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE)
# plotdata <- as.data.frame(plotdata)
# DateTimeStamp <- substr(filename, 1, 16)
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))


#### Showing and saving output ####
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
write.csv(plotdata, file = paste0(DateTimeStamp, "multispecies.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(mydatatotal, file = paste0(DateTimeStamp, "multispeciestotal.csv"),
          quote = FALSE, row.names = FALSE)

## Labels and limits for plots ##
labspecies <- paste(nspeciesset, "species")
names(labspecies) <- nspeciesset
labmodel <- c("Broken stick model", "Dominance preemption model")
names(labmodel) <- c(1, 2)
mylabeller <- labeller(nspecies = labspecies, modelcode = labmodel,
                       .default = label_both)

plotdata <- as.data.frame(plotdata)
mydatatotal <- as.data.frame(mydatatotal)
limitsfraction <- c(0, 1)
limitsgrowthrate <- c(min(plotdata[, "mingrowthrate"]),
                      max(plotdata[, "maxgrowthrate"]))
# Round the limits to one decimal place, while ensuring that all the data is
# within the rounded limits. 
limitsgrowthrate <- sign(limitsgrowthrate)*ceiling(abs(limitsgrowthrate)*10)/10
limitsgrowthratebinned <- sort(c(limitsgrowthrate, limitsgrowthrate/2, 0))


## Plot summary data for the calculated growth rates 
CreatePlot(fillvar = "mingrowthrate", filltitle = "Minimum growth rate",
           filltype = "binned", limits = limitsgrowthratebinned, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "minor")

CreatePlot(fillvar = "mingrowthrate", filltitle = "Minimum growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "minor")

CreatePlot(fillvar = "meangrowthrate", filltitle = "Mean growth rate",
           filltype = "binned", limits = limitsgrowthratebinned, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "minor")

CreatePlot(fillvar = "meangrowthrate", filltitle = "Mean growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "minor")

CreatePlot(fillvar = "maxgrowthrate", filltitle = "Max growth rate",
           filltype = "binned", limits = limitsgrowthratebinned, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "minor")

CreatePlot(fillvar = "maxgrowthrate", filltitle = "Max growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "minor")

selectmydatatotal <- filter(mydatatotal, iter == niter)
# Show how interactions affect the growth rates of the individual species
# required to obtain equilibrium. NOTE: costs do not affect growth rate, so
# differences between costs reflect stochastic effects.
# NOTE 2: data for all iterations is plotted on top of each other, so if all
# data is used, only the data for last iteration is visible.
# Using fillvar = "mean(growthrate)" does not work.
CreatePlot(dataplot = selectmydatatotal, fillvar = "growthrate",
           filltitle = "Growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "species + nspecies", facetx = "modelcode + cost + conjrate",
           diagional = "minor")


## Plot equilibrium characteristics for model without plasmids
CreatePlot(fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
CreatePlot(fillvar = "fracinfgrowth", filltitle = "Fraction infinite\ngrowth",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
CreatePlot(fillvar = "fraceqreached", filltitle = "Fraction equilibrium\nreached",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

## Plot total abundances of plasmid-free populations after perturbations for
# models without plasmids. Only abundances where perturbation did NOT lead to
# infinite growth are considered.
filltitle <- "Minimum total abundance of\nplasmid-free bacteria\nafter perturbation"
CreatePlot(fillvar = "minR", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
filltitle <- "Mean total abundance of\nplasmid-free bacteria\nafter perturbation"
CreatePlot(fillvar = "meanR", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
lengthx <- length(which(!is.na(plotdata[, "meanR"])))
plot(x = 1:lengthx, y = sort(log10(plotdata[, "meanR"])), type = "p", lwd = 2); grid()

ggplot(plotdata, aes(log10(meanR))) +
  geom_density(aes(fill = factor(nspecies)), show.legend = TRUE) +
  facet_grid(nspecies ~ modelcode)

filltitle <- "Median total abundance of\nplasmid-free bacteria\nafter perturbation"
CreatePlot(fillvar = "medianR", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
lengthx <- length(which(!is.na(plotdata[, "medianR"])))
plot(x = 1:lengthx, y = sort(log10(plotdata[, "medianR"]))); grid()

filltitle <- "Maximum total abundance of\nplasmid-free bacteria\nafter perturbation"
CreatePlot(fillvar = "maxR", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

## Plot equilibrium characteristics for model with plasmids
CreatePlot(fillvar = "fracstableconj",
           filltitle = "Fraction stable\nwith conjugation",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
CreatePlot(fillvar = "fracinfgrowthconj",
           filltitle = "Fraction infinite growth\nwith conjugation",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
CreatePlot(fillvar = "fraceqreachedconj",
           filltitle = "Fraction equilibrium\nreached with\nconjugation",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

## Plot total abundances of plasmid-free populations after perturbations for
# models with plasmids. Only abundances where perturbation did NOT lead to
# infinite growth are considered.
filltitle <- "Minimum total abundance of\nplasmid-free bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "minRconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
filltitle <- "Mean total abundance of\nplasmid-free bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "meanRconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

ggplot(plotdata, aes(log10(meanRconj))) +
  geom_density(aes(fill = factor(nspecies)), show.legend = TRUE) +
  facet_grid(nspecies + conjrate ~ modelcode + cost)

filltitle <- "Median total abundance of\nplasmid-free bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "medianRconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
filltitle <- "Maximum total abundance of\nplasmid-free bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "maxRconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")


## Plot total abundances of plasmid-bearing populations after perturbations for
# models with plasmids. Only abundances where perturbation did NOT lead to
# infinite growth are considered.
filltitle <- "Minimum total abundance of\nplasmid-bearing bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "minPconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
filltitle <- "Mean total abundance of\nplasmid-bearing bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "meanPconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
filltitle <- "Median total abundance of\nplasmid-bearing bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "medianPconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")
filltitle <- "Maximum total abundance of\nplasmid-bearing bacteria after\nperturbation with plasmids"
CreatePlot(fillvar = "maxPconj", filltitle = filltitle,
           filltype = "continuous", limits = NULL, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

## Plot other equilibrium characteristics for models with and without plasmids
CreatePlot(fillvar = "fracreal", filltitle = "Fraction real",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

CreatePlot(fillvar = "fracrealconj",
           filltitle = "Fraction real\nwith conjugation",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

CreatePlot(fillvar = "fracrep", filltitle = "Fraction repeated eigenvalues",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

CreatePlot(fillvar = "fracrepconj",
           filltitle = "Fraction repeated eigenvalues\nwith conjugation",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies + conjrate", facetx = "modelcode + cost",
           diagional = "both")

### To test plots without using CreatePlot() ###
# ggplot(data = plotdata, aes(x = intmean, y = selfintmean, fill = fracstable)) +
#   geom_raster() +
#   theme_bw(base_size = 15) +
#   scale_x_continuous() +
#   scale_y_continuous() +
#   scale_fill_viridis_c("Fraction stable", limits = limitsfraction) +
#   geom_abline(intercept = 0, slope = -1, col = "white", size = 1.1) +
#   geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1) +
#   coord_fixed(ratio = 1, expand = FALSE) +
#   theme(legend.position = "bottom") +
#   labs(x = "Mean interaction coefficient",
#        y = "Mean selfinteraction coefficient",
#        caption = paste(niter, "iterations")) +
#   facet_grid(nspecies ~ modelcode, labeller = mylabeller)


### Show species-specific relation of growth rate, intmean and selfintmean.
mydata <- as.data.frame(mydata)

# Check relation between abundance and growth rate for the different species.
nrow <- dim(mydatatotal)[1]
subsetmydatatotal <- filter(mydatatotal, near(iter, 1))

# Larger selfintmean leads to smaller growth rate, but effect becomes smaller
# when species are less abundant.
# Larger intmean leads to lower growth rate, effect becomes larger when species
# are less abundant
ggplot(data = subsetmydatatotal, aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = selfintmean), size = 1) +
  facet_grid(species + nspecies ~ modelcode + cost + conjrate, labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations"))
ggsave(paste0(DateTimeStamp, "growthrate1perspecies.png"))

ggplot(data = subsetmydatatotal, aes(x = selfintmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = intmean), size = 1) +
  facet_grid(species + nspecies ~ modelcode + cost + conjrate, labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations"))
ggsave(paste0(DateTimeStamp, "growthrate2perspecies.png"))


#### Comparing abundance models ####
comparingabuntotal <- NULL

for(nspecies in nspeciesset) {
  comparingabun <- data.frame(
    nspecies = as.factor(rep(nspecies, 2*nspecies)),
    species = as.factor(rep(1:nspecies, 2)),
    abun = c(brokenstick(nspecies = nspecies, totalabun = totalabun,
                         takelimit = TRUE),
             dompreempt(nspecies = nspecies, totalabun = totalabun,
                        takelimit = TRUE)),
    model = rep(c("brokenstick", "dompreempt"), each = nspecies))
  comparingabuntotal <- rbind(comparingabuntotal, comparingabun)
  
  plotabun <- ggplot(data = comparingabun, aes(x = species, y = abun, color = model)) +
    theme_bw() +
    scale_x_discrete(limits = factor(1:max(nspeciesset))) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "bottom") +
    labs(title = "Comparing abundance models: linear scale",
         x = "Species rank", y = "Species abundance") +
    geom_line(aes(group = model), size = 1.25) +
    geom_point(size = 2)
  print(plotabun)
  if(saveplots == TRUE) {
    filename <- paste0(DateTimeStamp, "compareabun", nspecies, "species.png")
    ggsave(filename)
  }
  
  plotabunlog <- ggplot(data = comparingabun, aes(x = species, y = abun, color = model)) +
    theme_bw() +
    scale_x_discrete(limits = factor(1:max(nspeciesset))) +
    scale_y_continuous(limits = c(NA, 1), trans = "log10") +
    theme(legend.position = "bottom") +
    labs(title = "Comparing abundance models: logarithmic scale",
         x = "Species rank", y = "Species abundance") +
    geom_line(aes(group = model), size = 1.25) +
    geom_point(size = 2)
  print(plotabunlog)
  if(saveplots == TRUE) {
    filename <- paste0(DateTimeStamp, "compareabun", nspecies, "specieslog.png")
    ggsave(filename)
  }
}
comparingabuntotal <- comparingabuntotal
comparingabuntotal[, "group"] <- paste0(comparingabuntotal[, "nspecies"],
                                             " species, ", comparingabuntotal[, "model"])

plotabun <- ggplot(data = comparingabuntotal, aes(x = species, y = abun, color = nspecies, lty = model)) +
  theme_bw(base_size = 14) +
  scale_x_discrete(limits = factor(1:max(nspeciesset))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "bottom") +
  labs(title = "Comparing abundance models: linear scale",
       x = "Species rank", y = "Species abundance") +
  geom_line(aes(group = group), size = 1.25) +
  geom_point(size = 2)
print(plotabun)
if(saveplots == TRUE) {
  filename <- paste0(DateTimeStamp, "compareabuntotal.png")
  ggsave(filename)
}

plotabunlog <- ggplot(data = comparingabuntotal, aes(x = species, y = abun, color = nspecies, lty = model)) +
  theme_bw(base_size = 14) +
  scale_x_discrete(limits = factor(1:max(nspeciesset))) +
  scale_y_continuous(limits = c(NA, 1), trans = "log10") +
  theme(legend.position = "bottom") +
  labs(title = "Comparing abundance models: logarithmic scale",
       x = "Species rank", y = "Species abundance") +
  geom_line(aes(group = group), size = 1.25) +
  geom_point(size = 2)
print(plotabunlog)
if(saveplots == TRUE) {
  filename <- paste0(DateTimeStamp, "compareabuntotallog.png")
  ggsave(filename)
}