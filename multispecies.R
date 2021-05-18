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
# Test if using larger timesteps speeds up simulations without affecting the
# result, since variable time-step method is used and most large changes occur
# early. For example, tstep <- 1;
# times <- c(0:300, seq(from = 300 + tstep, to = tmax, by = tstep))


#### Optionally to do ####

## General ##
# Formatting of text in messages, warnings, and graphs can be coded nicer using
# the stringr package.

# Some variables are named x -> xconj (maxabunR, maxabunRconj), others are named
# fraceigvalRep -> fraceigvalconjRep. That could be made more consistent.

# I could move the 'niter' argument to the top functions, such that for 1000
# iterations, rnorm is called only once to generate 1000*nspecies growthrates,
# instead of being called 1000 times to generate nspecies growthrates.
# Or use replicate(...) from the apply-family.

# To check if enough iterations are used: run several times for niter iterations,
# if the variation in fraction of stable equilibria is too large, use more
# iterations.

# Now I use nested loops, instead I could first create a tibble using
# tidyr::expand_grid(), and then use (l)/(m)apply / purrr:(p)map to iterate over
# all rows?

# The .groups argument in dplyr::summarise is experimental, so maybe should be
# dropped. However, excluding it leads to messages every time it is called.

# The matrix 'data' is converted to a tibble to efficiently get summary
# statistics. The tibble with summary statistics is then converted to a matrix
# to fill in part of the matrix plotdata. Using a tibble for 'data' from the
# start prevents some of this type-conversion, might make naming columns on
# assignment clearer, and enable use of the more efficient dplyr::bind_cols()
# instead of base::cbind(). However, originally using a data.frame instead of
# matrix to store data made progress really slow, so I should check how timing
# is affected when 'data', 'plotdata' or both are tibbles from the start.

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
# See May 2005 on other options for intraspecies interactions (all < 0, all -1, ...).
# Search literature for abundance models applied to a microbiome.

# Add logistic interaction in addition the the linear interaction currently
# implemented (see https://github.com/EgilFischer/FlockMicrobiome for code).

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

## createplot() ##
# The diagonals are not displayed correctly if the plotting area is non-square


#### Loading required libraries ####
library(deSolve)   # checkequilibrium and perturbequilibrium call ode() if
# showplot == TRUE and simulateinvasion == TRUE, respectively
library(dplyr)     # across(), group_by(), near(), summarise()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()
# On the pipe operator (%>%), see ?'%>%' and https://r4ds.had.co.nz/pipes.html


#### Settings and defining parameterspace ####

# Simulation settings
niter <- 10
niterintmat <- 1
simulateinvasion <- TRUE # If TRUE, simulations over time are performed
saveplots <- TRUE
smallstate <- 1e-20 # States are set to 0 if they become smaller than smallstate
smallchange <- 1e-10 # If the sum of absolute rates of change is equal to
# smallchange, equilibrium is assumed to be reached and integration is terminated

# Define parameter space
totalabun <- 1
nspeciesset <- c(2, 4, 6)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("brokenstick", "dompreempt")
intmeanset <- seq(from = -0.8, to = 0.6, by = 0.1)
selfintmeanset <- seq(from = -0.8, to = -0.3, by = 0.1)
costset <- c(0.01, 0.20)
conjrateset <- list(rep(0.01, maxnspecies), rep(0.05, maxnspecies),
                    rep(0.1, maxnspecies))
conjmattype <- "diffTax"  # NOTE: use of multiple values is NOT supported (yet)
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")

# Settings for testing code
niter <- 2
niterintmat <- 1
simulateinvasion <- TRUE # If TRUE, simulations over time are performed
saveplots <- TRUE
smallstate <- 1e-20
smallchange <- 1e-10
totalabun <- 1
nspeciesset <- c(2, 4, 6)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("brokenstick", "dompreempt")
intmeanset <- seq(from = -0.6, to = 0.6, by = 0.6)
selfintmeanset <- seq(from = -0.8, to = -0.5, by = 0.3)
costset <- c(0.01, 0.20)
conjrateset <- list(rep(0.01, maxnspecies), rep(0.1, maxnspecies))
conjmattype <- "diffTax"  # NOTE: use of multiple values is NOT supported (yet)
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")

# A matrix giving the taxonomic relatedness between the donor species (columns)
# and recipient species (rows). Although the matrix is symmetric, this ordering
# is needed to calculate the conjugation rates correctly. Only the submatrices
# taxmat[1:nspecies, 1:nspecies] are used if nspecies < max(nspeciesset).
# An example with E. coli, Klebsiella, two Erwinia species, and two Xanthomonas species:
taxmat <- matrix(c("SameSpecies", "SameFamily",  "SameOrder",   "SameOrder",   "SameClass",   "SameClass",
                   "SameFamily",  "SameSpecies", "SameOrder",   "SameOrder",   "SameClass",   "SameClass",
                   "SameOrder",   "SameOrder",   "SameSpecies", "SameSpecies", "SameClass",   "SameClass",
                   "SameOrder",   "SameOrder",   "SameSpecies", "SameSpecies", "SameClass",   "SameClass",
                   "SameClass",   "SameClass",   "SameClass",   "SameClass",   "SameSpecies", "SameSpecies",
                   "SameClass",   "SameClass",   "SameClass",   "SameClass",   "SameSpecies", "SameSpecies"),
                  nrow = 6, ncol = 6, byrow = TRUE)
taxmat <- matrix(rep("SameSpecies", maxnspecies^2),
                 nrow = maxnspecies, ncol = maxnspecies, byrow = TRUE)


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
# The roots should terminate the simulation when equilibrium is reached, but if
# one starts in an equilibrium, this will not happen.
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
               parms = list(growthrate = growthrate, intmat = intmat),
               rootfun = rootfun,
               events = list(func = eventfun, root = TRUE, terminalroot = c(1, 2)))
    ylim <- c(0, 1.1*max(out[, -1]))
    matplot.deSolve(out, ylim = ylim, lwd = 2,
                    lty = 1, ylab = "Abundance")
    grid()
  }
  return(atequilibrium)
}

# Get a matrix with conjugation rates, for different scenarios.
# Input:
#   - nspecies: integer indicating the number of species
#   - type: character vector determining if and how different conjugation rates
#     should be generated. Should be one of (1) 'diffDiag' to set all
#     intraspecies conjugation rates equal to the first element in conjrate, and
#     all off-diagonal elements equal to the second element in conjrate, or (2)
#     'diffTax' to set conjugation rates are based on taxonomic relatedness as
#     specified in taxmat.
#   - conjrate: numeric vector of length 2 giving intraspecies and interspecies
#     conjugation rates (if type = 'diffDiag'), or of length nspecies giving
#     the intraspecies conjugation rates for each species (if type = 'diffTax').
#   - taxmat: matrix where element in column n and row r gives taxonomic
#     relatedness between species n and r, as a character string (possible are
#     "SameSpecies", "SameFamily", "SameOrder", "SameClass", and "OtherClass").
#     This matrix is used to calculate the interspecies conjugation rates based
#     on the intraspecies conjugation rates provided in conjrate, and the
#     taxonomic relatedness as provided in taxmat. By default R fills matrices
#     by column, so byrow = TRUE should be used to obtain the same matrix as
#     the 'folded' vector. For example for 2 species:
#     taxmat <- matrix(c("SameSpecies", "SameFamily",
#                          "SameFamily", "SameSpecies"),
#                        nrow = nspecies, ncol = nspecies, byrow = TRUE)
# Return:
#   - conjmat: matrix where the element in column n and row r gives the
#     conjugation rate from species n to species r.
# Notes:
# - To obtain a matrix where all values are identical (i.e., where conjugation
#   is not dependent on relatedness), use type = 'diffDiag' with two identical
#   values for conjrate.
getconjmat <- function(nspecies, type, conjrate, taxmat = NULL) {
  switch(type,
         diffDiag = {
           stopifnot(near(length(conjrate), 2))
           # Fill matrix with interspecies conjugation rates
           conjmat <- matrix(rep(conjrate[2], nspecies^2),
                             nrow = nspecies, ncol = nspecies)
           # Put intraspecies conjugation rates on the diagonal
           diag(conjmat) <- conjrate[1]
         },
         diffTax = {
           stopifnot(near(length(conjrate), nspecies),
                     !is.null(taxmat), dim(taxmat)[1] == nspecies,
                     all(diag(taxmat) == "SameSpecies"),
                     isSymmetric.matrix(unname(taxmat)))
           # To obtain interspecies conjugation rates for the different levels
           # of taxonomic relatedness between donor and recipients, the
           # intraspecies conjugation rates are multiplied with the following
           # conversion factors.
           rateconv <- c(SameSpecies = 1, SameFamily = 2.7, SameOrder = 0.1,
                         SameClass = 0.05, OtherClass = 0.001)
           convmat <- matrix(NA, nrow = nspecies, ncol = nspecies)
           for(taxlevel in names(rateconv)) {
             convmat[which(taxmat == taxlevel)] <- rateconv[taxlevel]
           }
           # Multiply column n of convmat giving the conversion factors for
           # conjugation from donor species n to the different recipient species,
           # with element n of conjrate. Using t(t(convmat) * conjrate) is
           # faster for nspecies > 15
           conjmat <- convmat %*% diag(conjrate)
           conjmat
         },
         {
           warning("'type' should be 'diffDiag' or 'diffTax'.")
           conjmat <- NULL
         }
  )
  return(conjmat)
}

# Get equilibrium characteristics
geteqinfo <- function(model, abundance, intmat, growthrate,
                      cost = NULL, conjmat = NULL) {
  if(model == "gLV") {
    eigval <- eigen(
      x = jacobian.full(y = abundance, func = gLV,
                        parms = list(growthrate = growthrate, intmat = intmat)),
      symmetric = FALSE, only.values = TRUE)$values
  } else {
    eigval <- eigen(
      x = jacobian.full(y = abundance, func = gLVConj,
                        parms = list(growthrate = growthrate, intmat = intmat,
                                     cost = cost, conjmat = conjmat)),
      symmetric = FALSE, only.values = TRUE)$values
  }
  
  eigvalRep <- any(duplicated(eigval))
  # Using sort(eigval) to get eigenvalue with largest real part, because
  # max(eigval) does not work if eigval is complex. Complex values are sorted
  # first by the real part, then the imaginary part.
  eigval <- sort(eigval, decreasing = TRUE)[1]
  eigvalRe <- Re(eigval)
  eigvalIm <- Im(eigval)
  eigvalReSign <- sign(eigvalRe)
  eigvalImSign <- sign(eigvalIm)
  eqinfo <- c(eigvalRe, eigvalIm, eigvalReSign, eigvalImSign, eigvalRep)
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
                               tmax = 1e4, tstep = 0.1, showplot = TRUE,
                               verbose = FALSE, suppresswarninfgrowth = FALSE) {
  
  # Name abundances, set line type and colors, get derivatives of initial state
  # in the plasmid-free model.
  if(model == "gLV") {
    nspecies <- length(abundance)
    names(abundance) <- paste0(rep("R", nspecies), 1:nspecies)
    lty <- 1
    col <- mycol[1:nspecies]
    derivatives <- unlist(
      gLV(t = 0, n = abundance, parms = list(growthrate = growthrate, intmat = intmat))
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
      gLV(t = 0, n = abundance[1:nspecies],
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
  final <- tail(out, 1)
  timefinal <- final[, 1]
  abunfinaltemp <- final[, -1]
  names(abunfinaltemp) <- names(abundance)
  
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
                 paste(names(abunfinaltemp), "=", signif(abunfinaltemp), collapse = ", "))
        )
      }
    }
  }
  
  if(eqreached == 0 && infgrowth != 1) {
    warning("Equilibrium not reached. Increase tmax to prevent this? Final abundances were\n",
            paste(names(abunfinaltemp), "=", signif(abunfinaltemp), collapse = ", "))
  }
  
  if(showplot == TRUE) {
    subtitle <- paste0(abunmodel, ", intmean=", intmean, ", selfintmean=", selfintmean, ", cost=",
                       cost, ", conjratecode=", conjratecode)
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
    
    abschange <- abunfinaltemp - abuninit
    relchange <- abunfinaltemp / abuninit
    print("Species abundances, and their changes:", quote = FALSE)
    print(rbind(tested = abuninit, at_perturbation = abunpert, final = abunfinaltemp,
                absolute_change = abschange, relative_change = relchange))
  }
  
  if(model == "gLV") {
    abunfinal <- list(abunfinalR = abunfinaltemp,
                      abunfinalP = NULL,
                      infgrowth = infgrowth, eqreached = eqreached,
                      timefinal = timefinal)
  } else {
    abunfinal <- list(abunfinalR = abunfinaltemp[1:nspecies],
                      abunfinalP = abunfinaltemp[(nspecies + 1):(2*nspecies)],
                      infgrowth = infgrowth, eqreached = eqreached,
                      timefinal = timefinal)
  }
  return(abunfinal)
}

# 'Functions' to get summary statistics for the input data.
# Actually not functions, but lists of functions. They are called with 
# dplyr::summarise. I set na.rm = TRUE to get sensible values if not all values
# are NA, otherwise Inf, NaN, NA, and -Inf are returned.
getsummary4 <- list(
  min = ~min(.x, na.rm = TRUE),
  mean = ~mean(.x, na.rm = TRUE),
  median = ~median(.x, na.rm = TRUE),
  max = ~max(.x, na.rm = TRUE)
)
getfracnotzero <- list(frac = ~mean(.x != 0))

# Function to create plots
CreatePlot <- function(dataplot = plotdata, xvar = "intmean", yvar = "selfintmean",
                       fillvar, filltitle, filltype = "discrete", limits = NULL, 
                       title = NULL, subtitle = NULL,
                       labx = "Mean interaction coefficient",
                       laby = "Mean selfinteraction coefficient",
                       tag = NULL, addstamp = FALSE, diagonal = "none",
                       facetx = "abunmodelcode + cost",
                       facety = "nspecies + conjratecode",
                       as.table = TRUE,
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
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar)) + 
    theme_bw(base_size = base_size) +
    geom_raster() +
    scale_x_continuous() +
    scale_y_continuous() +
    coord_fixed(ratio = 1, expand = FALSE) +
    facet_grid(as.formula(paste(facety, "~", facetx)), as.table = as.table,
               labeller = mylabeller) +
    theme(legend.position = "bottom") +
    labs(title = title, subtitle = subtitle,
         x = labx, y = laby, caption = caption, tag = tag)
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
  if(diagonal == "both" || diagonal == "major") {
    p <- p + geom_abline(intercept = 0, slope = -1, col = "white", size = 1.1)
  }
  if(diagonal == "both" || diagonal == "minor") {
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
# geteqinfo(model = "gLV", abundance = abunbrokenstick, intmat = intmat1, growthrate = growthratebrokenstick)
# geteqinfo(model = "gLV", abundance = abundompreempt, intmat = intmat1, growthrate = growthratedompreempt)


#### Running the simulations ####
set.seed(seed = 314, kind = "default", normal.kind = "default", sample.kind = "default")

# Create matrix to store data
nrowplotdata <- prod(lengths(list(nspeciesset, abunmodelset, intmeanset,
                                  selfintmeanset, costset, conjrateset),
                             use.names = FALSE))
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata,
                   ncol = if(simulateinvasion == TRUE) {47 + 12*maxnspecies} else {25})
nrowdatatotal <- prod(lengths(list(abunmodelset,intmeanset, selfintmeanset,
                                   costset, conjrateset), use.names = FALSE))*niter*sum(nspeciesset)
datatotal <- matrix(data = NA, nrow = nrowdatatotal, ncol = 33 + 3*maxnspecies)
indexdatatotal <- 1

# Run simulations
rowindexplotdata <- 1
rowindexdata <- 1
for(nspecies in nspeciesset) {
  
  for(abunmodel in abunmodelset) {
    if(abunmodel == "brokenstick") {
      abundance <- brokenstick(nspecies = nspecies, totalabun = totalabun,
                               takelimit = TRUE)
      # Using a number instead of a name, to prevent type-conversion when
      # storing it in matrices.
      abunmodelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt(nspecies = nspecies, totalabun = totalabun,
                              takelimit = TRUE)
      abunmodelcode <- 2
    }
    
    for(intmean in intmeanset) {
      print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
                   ", intmean = ", intmean,
                   ": started at ", Sys.time()), quote = FALSE)
      
      for(selfintmean in selfintmeanset) {
        # print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
        #              ", intmean = ", intmean, ", selfintmean = ", selfintmean,
        #              ": started at ", Sys.time()), quote = FALSE)
        nrowdata <- niter * nspecies * length(costset) * length(conjrateset)
        data <- matrix(data = NA, nrow = nrowdata, ncol = 33 + 3*maxnspecies)
        indexdata <- 1
        abunR <- rep(NA, maxnspecies)
        abunRconj <- rep(NA, maxnspecies)
        abunPconj <- rep(NA, maxnspecies)
        
        for(iter in 1:niter) {
          stableeq <- FALSE
          iterintmat <- 0
          
          # Create a combination of interaction matrix and growth rate that
          # results in a stable plasmid-free equilibrium in the model
          # without conjugation
          while(stableeq == FALSE && iterintmat < niterintmat) {
            intmat <- getintmat(nspecies = nspecies,
                                intmean = intmean, selfintmean = selfintmean)
            growthrate <- getgrowthrate(abundance = abundance, intmat = intmat)
            eqinfo <- geteqinfo(model = "gLV", abundance = abundance,
                                intmat = intmat, growthrate = growthrate)
            if(eqinfo[1] < 0) {
              stableeq <- TRUE
            }
            iterintmat <- iterintmat + 1
          }
          if(stableeq == FALSE) {
            warning(paste("No stable equilibrium has been found in", niterintmat, "attempts."))
          }
          
          # Simulate invasion with plasmid-free bacteria of plasmid-free
          # equilibrium in model without conjugation (i.e., test internal
          # stability)
          if(simulateinvasion == TRUE) {
            if(stableeq == FALSE) {
              abunfinal <- perturbequilibrium(abundance = abundance, intmat = intmat,
                                              growthrate = growthrate, cost = cost,
                                              conjmat = conjmat,
                                              model = "gLV", pertpop = "all", tmax = 1e4,
                                              showplot = FALSE, verbose = FALSE,
                                              suppresswarninfgrowth = TRUE)
            } else {
              # No need for simulations if equilibrium is stable
              abunfinal <- list(abunfinalR = abundance,
                                abunfinalP = NULL,
                                infgrowth = 0, eqreached = 1,
                                timefinal = 1)
            }
            infgrowth <- abunfinal$infgrowth
            eqreached <- abunfinal$eqreached
            timefinal <- abunfinal$timefinal
            
            # It does not make sense to store abundances in case of infinite
            # growth or if equilibrium is not reached, so record those as NA
            if(eqreached == 1) {
              abunR[1:nspecies] <- abunfinal$abunfinalR
              abunRtotal <- sum(abunfinal$abunfinalR)
            } else {
              abunR[1:nspecies] <- NA
              abunRtotal <- NA
            }
          } else {
            # No simulations over time performed, so set values to NA
            infgrowth <- NA
            eqreached <- NA
            timefinal <- NA
            abunR[1:nspecies] <- NA
            abunRtotal <- NA
          }
          
          for(cost in costset) {
            conjratecode <- 0
            
            for(conjrate in conjrateset) {
              conjratecode <- conjratecode + 1
              conjratensp <- conjrate[1:nspecies]
              taxmatnsp <- taxmat[1:nspecies, 1:nspecies]
              conjmat <- getconjmat(nspecies = nspecies, type = conjmattype,
                                    conjrate = conjratensp,
                                    taxmat = taxmatnsp)
              # Get equilibrium characteristics for plasmid-free equilibrium in
              # the model with conjugation
              eqinfoconj <- geteqinfo(model = "gLVconj",
                                      abundance = c(abundance, rep(0, nspecies)),
                                      intmat = intmat, growthrate = growthrate,
                                      cost = cost, conjmat = conjmat)
              
              # To simulate invasion of the plasmid-free equilibrium with plasmids,
              # the abundances of the plasmid-free populations have to be appended
              # to the abundances of the plasmid-bearing populations
              if(simulateinvasion == TRUE) {
                if(eqinfoconj[1] >= 0) {
                  abunfinalconj <- perturbequilibrium(abundance = c(abundance, rep(0, nspecies)),
                                                      intmat = intmat, growthrate = growthrate,
                                                      cost = cost, conjmat = conjmat,
                                                      model = "gLVConj", pertpop = "P1", tmax = 1e4,
                                                      showplot = FALSE, verbose = FALSE,
                                                      suppresswarninfgrowth = TRUE)
                } else {
                  # No need for simulations if equilibrium is stable
                  abunfinal <- list(abunfinalR = abundance,
                                    abunfinalP = rep(0, nspecies),
                                    infgrowth = 0, eqreached = 1,
                                    timefinal = 1)
                }
                infgrowthconj <- abunfinalconj$infgrowth
                eqreachedconj <- abunfinalconj$eqreached
                timefinalconj <- abunfinalconj$timefinal
                
                # It does not make sense to store abundances in case of infinite
                # growth or if equilibrium is not reached, so record those as NA
                if(eqreachedconj == 1) {
                  abunRconj[1:nspecies] <- abunfinalconj$abunfinalR
                  abunPconj[1:nspecies] <- abunfinalconj$abunfinalP
                  abunRtotalconj <- sum(abunfinalconj$abunfinalR)
                  abunPtotalconj <- sum(abunfinalconj$abunfinalP)
                } else {
                  abunRconj[1:nspecies] <- NA
                  abunPconj[1:nspecies] <- NA
                  abunRtotalconj <- NA
                  abunPtotalconj <- NA
                }
              } else {
                # No simulations over time performed, so set values to NA
                infgrowthconj <- NA
                eqreachedconj <- NA
                timefinalconj <- NA
                abunRconj[1:nspecies] <- NA
                abunPconj[1:nspecies] <- NA
                abunRtotalconj <- NA
                abunPtotalconj <- NA
              }
              indexdatanew <- indexdata + nspecies
              
              data[indexdata:(indexdatanew - 1), ] <- cbind(
                niter, nspecies, abunmodelcode, intmean, selfintmean,
                cost, conjratecode, iter, 1:nspecies, abundance,
                diag(intmat), c(growthrate), iterintmat,
                matrix(rep(eqinfo, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(eqinfoconj, nspecies), nrow = nspecies, byrow = TRUE),
                infgrowth, infgrowthconj, eqreached, eqreachedconj,
                timefinal, timefinalconj, abunRtotal, abunRtotalconj, abunPtotalconj,
                abunRtotalconj/(abunRtotalconj + abunPtotalconj),
                matrix(rep(abunR, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(abunRconj, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(abunPconj, nspecies), nrow = nspecies, byrow = TRUE))
              indexdata <- indexdatanew
            }
          }
        }
        datatotal[indexdatatotal:(indexdatatotal + nrowdata - 1), ] <- data
        indexdatatotal <- indexdatatotal + nrowdata
        
        colnames(data) <- c("niter", "nspecies", "abunmodelcode",
                            "intmean", "selfintmean", "cost", "conjratecode",
                            "iter", "species", "abundance",
                            "selfintdata", "growthrate",
                            "iterintmat",
                            "eigvalRe", "eigvalIm",
                            "eigvalReSign", "eigvalImSign", "eigvalRep",
                            "eigvalconjRe", "eigvalconjIm",
                            "eigvalconjReSign", "eigvalconjImSign", "eigvalconjRep",
                            "infgrowth", "infgrowthconj",
                            "eqreached", "eqreachedconj",
                            "timefinal", "timefinalconj",
                            "abunRtotal", "abunRtotalconj", "abunPtotalconj",
                            "fracRtotalconj",
                            paste0("abunRsp", 1:maxnspecies),
                            paste0("abunRconjsp", 1:maxnspecies),
                            paste0("abunPconjsp", 1:maxnspecies))
        
        # Get summary data which do not depend on simulated invasion for all
        # combinations of costs and conjugation rates
        summarydata <- as_tibble(data) %>%
          group_by(cost, conjratecode) %>%
          summarise(
            across(c(selfintdata, growthrate, iterintmat),
                   getsummary4, .names = "{.col}{.fn}"),
            fracstable = mean(eigvalRe < 0),
            fracstableconj = mean(eigvalconjRe < 0),
            fracreal = mean(eigvalIm == 0),
            fracrealconj = mean(eigvalconjIm == 0),
            across(c(eigvalRep, eigvalconjRep), getfracnotzero,
                   .names = "{.fn}{.col}"),
            .groups = "drop"
          )
        
        if(simulateinvasion == TRUE) {
          # If invasion was simulated, get summary of generated data for all
          # combinations of costs and conjugation rates
          summarydatasimulation <- as_tibble(data) %>%
            group_by(cost, conjratecode) %>%
            summarise(
              across(c(starts_with("abunR"), starts_with("abunP"), fracRtotalconj),
                     getsummary4, .names = "{.col}{.fn}"),
              across(c(infgrowth, infgrowthconj, eqreached, eqreachedconj),
                     getfracnotzero, .names = "{.col}{.fn}"),
              timefinalmedian = median(timefinal),
              timefinalconjmedian = median(timefinalconj),
              .groups = "drop"
            )
          summarydata <- full_join(summarydata, summarydatasimulation,
                                by = c("cost", "conjratecode"))
        }
        
        rowindexplotdatanew <- rowindexplotdata + length(costset) * length(conjrateset)
        plotdata[rowindexplotdata:(rowindexplotdatanew - 1), ] <- as.matrix.data.frame(
          tibble(niter, nspecies, abunmodelcode, intmean, selfintmean, summarydata))
        rowindexplotdata <- rowindexplotdatanew
      }
    }
  }
}
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)
colnames(plotdata) <- c("niter", "nspecies", "abunmodelcode",
                        "intmean", "selfintmean", colnames(summarydata))
colnames(datatotal) <- colnames(data)

#### Saving output to .csv files ####
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
write.csv(plotdata, file = paste0(DateTimeStamp, "multispecies.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(datatotal, file = paste0(DateTimeStamp, "multispeciestotal.csv"),
          quote = FALSE, row.names = FALSE)

# Saving settings
names(conjrateset) <- paste0("conjrateset", 1:length(conjrateset))
rownames(taxmat) <- paste0("recipientsp", 1:dim(taxmat)[1])
colnames(taxmat) <- paste0("donorsp", 1:dim(taxmat)[1])
settings <- c(list(niter = niter, niterintmat = niterintmat,
                   simulateinvasion = simulateinvasion,
                   smallstate = smallstate, smallchange = smallchange,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   costset = costset), conjrateset,
              list(conjmattype = conjmattype, taxmat = taxmat))
for(index in 1:length(settings)) {
  write.table(t(as.data.frame(settings[index])), 
              paste0(DateTimeStamp, "settings.csv"), append = TRUE,
              quote = FALSE, sep = ",", col.names = FALSE)
}


#### Reading previously saved data from a .csv-file ####
## To read data from csv-file, uncomment this section and fill in the 
# needed datetimestamp
# filename <- "2021_05_04_17_44multispecies.csv"
# plotdata <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE)
# If plotdata has only one column, probably a semicolon instead of a comma is
# used as separator in .csv-files. So read the file again.
# if(dim(plotdata)[2] == 1) {
#   plotdata <- read.csv(filename, header = TRUE, sep = ";", quote = "\"",
#                        dec = ".", stringsAsFactors = FALSE)
# }
# plotdata <- as.data.frame(plotdata)
# DateTimeStamp <- substr(filename, 1, 16)
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))


#### Labels and limits for plots ####
labspecies <- paste(nspeciesset, "species")
names(labspecies) <- nspeciesset
labmodel <- c("Broken stick model", "Dominance preemption model")
names(labmodel) <- c(1, 2)
mylabeller <- labeller(nspecies = labspecies, abunmodelcode = labmodel,
                       .default = label_both)

plotdata <- as.data.frame(plotdata)
datatotal <- as.data.frame(datatotal)
limitsfraction <- c(0, 1)
# Round the limits to one decimal place, while ensuring that all the data is
# within the rounded limits. 
limitsgrowthrate <- c(floor(min(plotdata[, "growthratemin"])*10)/10,
                      ceiling(max(plotdata[, "growthratemax"])*10)/10)


#### To test plots without using CreatePlot() ####
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
#   facet_grid(nspecies ~ abunmodelcode, labeller = mylabeller)

#### Plot output ####

## Compare equilibrium characteristics for the models without and with plasmids.
# If invasion has been simulated, data on infinite growth and reaching
# equilibrium after perturbation is also plotted.
CreatePlot(fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracstableconj",
           filltitle = "Fraction stable\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)

# Show dominant eigenvalues
datatotalfilteredspecies <- filter(datatotal, near(species, 1))
datatotalfilteredspeciesiter <- filter(datatotalfilteredspecies, near(iter, 1))
limitseigvalRe <- range(c(datatotalfilteredspeciesiter[, "eigvalRe"],
                       datatotalfilteredspeciesiter[, "eigvalconjRe"]))
limitseigvalIm <- range(c(datatotalfilteredspeciesiter[, "eigvalIm"],
                          datatotalfilteredspeciesiter[, "eigvalconjIm"]))
CreatePlot(dataplot = datatotalfilteredspeciesiter,
           fillvar = "eigvalRe",
           filltitle = "Real part of\ndominant eigenvalue",
           filltype = "continuous", limits = limitseigvalRe)
CreatePlot(dataplot = datatotalfilteredspeciesiter,
           fillvar = "eigvalconjRe",
           filltitle = "Real part of\ndominant eigenvalue\nwith conjugation",
           filltype = "continuous", limits = limitseigvalRe)

limitseigvalRe <- range(c(datatotalfilteredspecies[, "eigvalRe"],
                          datatotalfilteredspecies[, "eigvalconjRe"]))
limitseigvalIm <- range(c(datatotalfilteredspecies[, "eigvalIm"],
                          datatotalfilteredspecies[, "eigvalconjIm"]))
ggplot(data = datatotalfilteredspecies,
       aes(x = eigvalRe, y = eigvalIm,
           color = as.factor(eigvalReSign))) +
  scale_x_continuous(limits = limitseigvalRe) +
  scale_y_continuous(limits = limitseigvalIm) +
  geom_point() +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(x = "Real part dominant eigenvalue",
       y = "Imaginary part dominant eigenvalue",
       caption = paste(niter, "iterations")) +
  facet_grid(rows = vars(nspecies), cols = vars(abunmodelcode),
             labeller = mylabeller)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "eigenvaluesdistr.png"))
}

ggplot(data = datatotalfilteredspecies,
       aes(x = eigvalconjRe, y = eigvalconjIm,
           color = as.factor(eigvalconjReSign))) +
  scale_x_continuous(limits = limitseigvalRe) +
  scale_y_continuous(limits = limitseigvalIm) +
  geom_point() +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(x = "Real part dominant eigenvalue (with conjugation)",
       y = "Imaginary part dominant eigenvalue (with conj.)",
       caption = paste(niter, "iterations")) +
  facet_grid(rows = vars(nspecies), cols = vars(abunmodelcode),
             labeller = mylabeller)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "eigenvaluesdistrconj.png"))
}

if(simulateinvasion == TRUE) {
  CreatePlot(fillvar = "infgrowthfrac", filltitle = "Fraction infinite\ngrowth",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedfrac", filltitle = "Fraction equilibrium\nreached",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "infgrowthconjfrac",
             filltitle = "Fraction infinite growth\nwith conjugation",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedconjfrac",
             filltitle = "Fraction equilibrium\nreached with\nconjugation",
             filltype = "continuous", limits = limitsfraction)
}
CreatePlot(fillvar = "fracreal", filltitle = "Fraction real",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracrealconj",
           filltitle = "Fraction real\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fraceigvalRep", filltitle = "Fraction repeated eigenvalues",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fraceigvalconjRep",
           filltitle = "Fraction repeated eigenvalues\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)

## Growth rates
# Plot summary data for the calculated growth rates 
CreatePlot(fillvar = "growthratemin", filltitle = "Minimum growth rate",
           filltype = "continuous", limits = limitsgrowthrate)
CreatePlot(fillvar = "growthratemean", filltitle = "Mean growth rate",
           filltype = "continuous", limits = limitsgrowthrate)
CreatePlot(fillvar = "growthratemedian", filltitle = "Median growth rate",
           filltype = "continuous", limits = limitsgrowthrate)
CreatePlot(fillvar = "growthratemax", filltitle = "Max growth rate",
           filltype = "continuous", limits = limitsgrowthrate)

# Show the relation of interactions and species-specific growth rate required to
# obtain an equilibrium. Costs and conjugation rate do not affect growth rate,
# so data has been filtered to have only one value for them.
datatotalfiltercostconj <- filter(datatotal, near(cost, costset[1]),
                                  near(conjratecode, 1))

# Larger intmean leads to lower growth rate, effect becomes larger when species
# are less abundant
ggplot(data = datatotalfiltercostconj, aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = selfintmean), size = 1) +
  facet_grid(species + nspecies ~ abunmodelcode, labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations"))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthrate1perspecies.png"))
}

ggplot(data = datatotalfiltercostconj, aes(x = selfintmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = intmean), size = 1) +
  facet_grid(species + nspecies ~ abunmodelcode, labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations"))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthrate2perspecies.png"))
}


## Plot summary data on the number of iterations in creating intmat needed to
# find a stable equilibrium with the model without plasmids
CreatePlot(fillvar = "iterintmatmin", filltitle =
             paste("Minimum number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(0, niterintmat))
CreatePlot(fillvar = "iterintmatmean", filltitle = 
             paste("Mean number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(0, niterintmat))
CreatePlot(fillvar = "iterintmatmedian", filltitle = 
             paste("Median number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(0, niterintmat))
CreatePlot(fillvar = "iterintmatmax", filltitle = 
             paste("Maximum number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(0, niterintmat))

if(simulateinvasion == TRUE) {
  CreatePlot(fillvar = "timefinalmedian", filltitle = "Median time",
             filltype = "continuous", title = "Time after perturbation",
             subtitle = "Perturbation with plasmid-free bacteria")
  CreatePlot(fillvar = "timefinalconjmedian", filltitle = "Median time",
             filltype = "continuous", title = "Time after perturbation",
             subtitle = "Perturbation with plasmid-bearing bacteria")
  
  ## Plot of total abundances of plasmid-free populations after perturbations in
  # models without plasmids. Only abundances where equilibrium was reached are
  # considered. Although costs and conjugation rates do not influence the
  # outcome, they are included as facets to facilitate comparison with plots of
  # abundances after perturbation with plasmid-bearing bacteria.
  title <- "Total abundances after perturbation"
  subtitle <- "Perturbation with plasmid-free bacteria"
  CreatePlot(fillvar = "abunRtotalmin", filltitle = "Minimum of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunRtotalmean", filltitle = "Mean of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunRtotalmedian", filltitle = "Median of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunRtotalmax", filltitle = "Maximum of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)

  ## Plot total abundances of plasmid-free populations after perturbations in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  title <- "Total abundances after perturbation"
  subtitle <- "Perturbation with plasmid-bearing bacteria"
  CreatePlot(fillvar = "abunRtotalconjmin", filltitle = "Minimum of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunRtotalconjmean", filltitle = "Mean of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunRtotalconjmedian", filltitle = "Median of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunRtotalconjmax", filltitle = "Maximum of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)  
  
  ## Plot total abundances of plasmid-bearing populations after perturbations for
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  title <- "Total abundances after perturbation"
  subtitle <- "Perturbation with plasmid-bearing bacteria"
  CreatePlot(fillvar = "abunPtotalconjmin", filltitle = "Minimum of plasmid-\nbearing bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunPtotalconjmean", filltitle = "Mean of plasmid-\nbearing bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunPtotalconjmedian", filltitle = "Median of plasmid-\nbearing bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "abunPtotalconjmax", filltitle = "Maximum of plasmid-\nbearing bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)

  ## Plot fraction of plasmid-free bacteria after perturbations for
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  title <- "Fraction plasmid-free bacteria after perturbation"
  subtitle <- "Perturbation with plasmid-bearing bacteria"
  CreatePlot(fillvar = "fracRtotalconjmin",
             filltitle = "Minimum fraction of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "fracRtotalconjmean",
             filltitle = "Mean fraction of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "fracRtotalconjmedian",
             filltitle = "Median fraction of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
  CreatePlot(fillvar = "fracRtotalconjmax",
             filltitle = "Maximum fraction of plasmid-\nfree bacteria",
             filltype = "continuous", title = title, subtitle = subtitle)
}

## Compare abundance models ##
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
    scale_x_discrete(limits = factor(1:maxnspecies)) +
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
    scale_x_discrete(limits = factor(1:maxnspecies)) +
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
  scale_x_discrete(limits = factor(1:maxnspecies)) +
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
  scale_x_discrete(limits = factor(1:maxnspecies)) +
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
