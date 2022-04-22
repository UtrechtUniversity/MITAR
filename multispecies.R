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

# Roberts MG, Heesterbeek JAP. 2021. Infection dynamics in ecosystems: on the
# interaction between red and grey squirrels, pox virus, pine martens and trees.
# Journal of the Royal Society, Interface 18(183):20210551.

# Tokeshi M. 1990. Niche apportionment or random assortment: species abundance
# patterns revisited. Journal of animal ecology 59(3):1129-1146.

# Wickham H, Grolemund G. 2017. R for data science: import, tidy, transform,
# visualize, and model data. Online version: https://r4ds.had.co.nz/index.html


#### To do ####

## geteqinfo() ##
# Currently only the sign and complex part of the largest eigenvalue is used
# to determine the type of equilibrium. However, sometimes the largest
# eigenvalue does not have a complex part when (some of) the other eigenvalues
# do have a complex part. If this affects the equilibrium this should be taken
# into account.
# I check for repeated eigenvalues, but note that the pair a +/- bi are not
# repeated eigenvalues. See p. 133 of Edelstein-Keshet 2005 and section 10.4.3.3 from
# https://eng.libretexts.org/Bookshelves/Industrial_and_Systems_Engineering/Book%3A_Chemical_Process_Dynamics_and_Controls_(Woolf)


#### Optionally to do ####

## General ##
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

## getintmat() ##
# Add logistic interaction in addition the the linear interaction currently
# implemented (see https://github.com/EgilFischer/FlockMicrobiome for code).

## geteqinfo() ##
# Determining the sign of the real and complex parts can be done vectorised
# outside the loop

## perturbequilibrium() ##
# I could try if supplying the analytic Jacobian to the solver speeds up the
# integration. See Box 1 in Lischke 2017.

## CreatePlot() ##
# The diagonals are not displayed correctly if the plotting area is non-square

# Adjust CreatePlot() to be able to take a list of additional arguments, such
# that these arguments are included in the saved plots. Now as a workaround set
# save to FALSE and used ggsave() to include the added guides arguments.


#### Loading required packages ####
library(deSolve)   # checkequilibrium and perturbequilibrium call ode() if
# showplot == TRUE and simulateinvasion == TRUE, respectively
library(dplyr)     # across(), group_by(), near(), summarise()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()
# On the pipe operator (%>%), see ?'%>%' and Ch. 18 'Pipes' in Wickham 2017.


#### Settings and defining parameter space ####
# If simulateinvasion == TRUE, simulations over time are performed
# If states become smaller than smallstate during the integration, they are set
# to 0
# If the sum of absolute rates of change is equal to smallchange, equilibrium is
#   assumed to be reached and the integration is terminated
# See the functions that use the arguments for more detailed info


## Basis parameter set
bifurparms <- FALSE
saveplots <- TRUE
smallstate <- 1e-3
finalsmallstate <- 1
smallchange <- 1e-2
totalabun <- 1e11
nspeciesset <- c(2, 8, 16)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
costset <- c(0.05, 0.09)
costtype <- "absolute"
conjrateset <- list(rep(1e-12, maxnspecies))
# If taxmattype is SameSpecies, the conjugation rate is the same for all
# populations. If taxmattype is PInOtherClass, the interspecies conjugation rate
# from and to the initially plasmid-bearing population is reduced a 1000-fold
taxmattypeset <- c("SameSpecies", "PInOtherClass")
# Introduce plasmid in the most abundant species (PInMostAbun == TRUE) or in the
# least abundant species (PInMostAbun == FALSE)
PInMostAbun <- TRUE
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")

## Parameters for detailed local stability analysis, not simulating invasion
niter <- 100
niterintmat <- 1
simulateinvasion <- FALSE
intmeanset <- seq(from = -1e-11, to = 1e-11, by = 1e-12)
selfintmeanset <- seq(from = -1e-11, to = 0, by = 1e-12)

## Smaller parameter set to simulate invasion
niter <- 50
niterintmat <- 1
simulateinvasion <- TRUE
intmeanset <- seq(from = -1e-11, to = 1e-11, by = 1e-12)
selfintmeanset <- seq(from = -1e-11, to = 0, by = 1e-12)

## Testset
bifurparms <- FALSE
saveplots <- TRUE
smallstate <- 1e-3
finalsmallstate <- 1
smallchange <- 1e-2
totalabun <- 1e11
nspeciesset <- c(2, 8, 16)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
costset <- c(0.05, 0.09)
costtype <- "absolute"
conjrateset <- list(rep(1e-13, maxnspecies), rep(1e-12, maxnspecies))
# If taxmattype is SameSpecies, the conjugation rate is the same for all
# populations. If taxmattype is PInOtherClass, the interspecies conjugation rate
# from and to the initially plasmid-bearing population is reduced a 1000-fold
taxmattypeset <- c("SameSpecies", "PInOtherClass")
# Introduce plasmid in the most abundant species (PInMostAbun == TRUE) or in the
# least abundant species (PInMostAbun == FALSE)
PInMostAbun <- TRUE
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")
niter <- 10
niterintmat <- 1
simulateinvasion <- TRUE
intmeanset <- seq(from = -0.9e-11, to = 0.9e-11, by = 3e-12)
selfintmeanset <- seq(from = -0.9e-11, to = 0, by = 3e-12)

## Small parameter set to test code
bifurparms <- FALSE
niter <- 2
niterintmat <- 1
simulateinvasion <- TRUE
saveplots <- TRUE
smallstate <- 1e-3
finalsmallstate <- 1
smallchange <- 1e-2
totalabun <- 1e11
nspeciesset <- c(2, 16)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
intmeanset <- c(-6e-12, 6e-12)
selfintmeanset <- c(-1.2e-11, -6e-12)
costset <- c(0.05, 0.09)
costtype <- "absolute"
conjrateset <- list(rep(1e-13, maxnspecies), rep(1e-12, maxnspecies))
# If taxmattype is SameSpecies, the conjugation rate is the same for all
# populations. If taxmattype is PInOtherClass, the interspecies conjugation rate
# from and to the initially plasmid-bearing population is reduced a 1000-fold
taxmattypeset <- c("SameSpecies", "PInOtherClass")
# Introduce plasmid in the most abundant species (PInMostAbun == TRUE) or in the
# least abundant species (PInMostAbun == FALSE)
PInMostAbun <- TRUE
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")

## Parameter set to create bifurcation-like plots showing the border of
# epidemiological stability in the conjugation rate/cost space
bifurparms <- TRUE
saveplots <- TRUE
smallstate <- 1e-3
finalsmallstate <- 1
smallchange <- 1e-2
totalabun <- 1e11
nspeciesset <- c(2, 8, 16)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
costset <- seq(from = 0, to = 0.2, by = 0.002)
costtype <- "absolute"
seqconjrate <- 10^seq(from = -13, to = -11, by = 0.02)
conjrateset <- NULL
for(conjrate in seqconjrate) {
  conjrateset <- c(conjrateset, list(rep(conjrate, maxnspecies)))
}
# If taxmattype is SameSpecies, the conjugation rate is the same for all
# populations. If taxmattype is PInOtherClass, the interspecies conjugation rate
# from and to the initially plasmid-bearing population is reduced a 1000-fold
taxmattypeset <- c("SameSpecies", "PInOtherClass")
# Introduce plasmid in the most abundant species (PInMostAbun == TRUE) or in the
# least abundant species (PInMostAbun == FALSE)
PInMostAbun <- TRUE
mycol <- c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
           "darkorange", "green1", "yellow", "hotpink")
niter <- 1
niterintmat <- 1
simulateinvasion <- FALSE
intmeanset <- c(-1e-11, 0, 1e-11) # seq(from = -1e-11, to = 1e-11, by = 1e-11)
selfintmeanset <- -0.5e-11 # seq(from = -1e-11, to = 0, by = 1e-11)


# Note: to simulate that each species belongs to a different class, an
# additional taxmattype has to be added. Then only the interspecies conjugation
# rates should be reduced, as the intraspecies conjugation rates of all species
# are still those given in conjrateset. This unchanged intraspecies conjugation
# rate makes it different from reducing conjrateset 1000-fold and using
# 'taxmatsame' as taxmattype.


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
# In matrix-notation this becomes: dn/dt = r*n*(1 - (1/K) * intmat %*% n).
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
# (= MacArthur fraction model) in the description of Tokeshi (1990). With the
# default takelimits = TRUE, this is done niterabun times and the mean abundance
# for each species is returned.
brokenstick <- function(nspecies, totalabun, takelimit = TRUE, niterabun = 1000) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
            length(totalabun) == 1, totalabun > 0)
  if(takelimit == TRUE) {
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
# subsequent species occupies more than half of the remainder (Tokeshi 1990).
# Over many iterations, each species preempts on average (0.5 + 1)/2 = 0.75 of
# the remainder, such that the abundances converge to the geometric series where
# the abundance of species i is given by totalabun*k*(1 - k)^(i - 1) with k =
# 0.75. The geometric model is used instead of the dominance preemption model if
# takelimit = TRUE, which is the default. The abundances are divided by
# 1 - (1 - 0.75)^nspecies to obtain the user-defined total abundance (Tokeshi 1990).
dompreempt <- function(nspecies, totalabun, takelimit = TRUE) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
            length(totalabun) == 1, totalabun > 0)
  abun <- rep(NA, nspecies)
  if(takelimit == TRUE) {
    abun <- totalabun*0.75*0.25^((1:nspecies) - 1)
  } else {
    remainingabun <- totalabun
    for(speciesindex in 1:nspecies) {
      abun.temp <- runif(1, min = 0.5, max = 1)*remainingabun
      abun[speciesindex] <- abun.temp
      remainingabun <- remainingabun - abun.temp
    }
  }
  # Scaling abundances to get specified total abundance
  abun <- abun / (1 - (0.25^nspecies))
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
# maximum of the range, or provide the standard deviation of the normal
# distribution as zero. The other arguments specify the distributions from which
# interaction coefficients are drawn.
getintmat <- function(nspecies, sparsity = 0,
                      intdistr = "normal", intmean = 0, intsd = 5e-12,
                      intrange = c(-1.2e-11, 1.2e-11),
                      selfintdistr = "normal",
                      selfintmean = -6e-13, selfintsd = 9e-12,
                      selfintrange = c(-1.2e-11, 0)) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
            is.numeric(sparsity), length(sparsity) == 1,
            sparsity >= 0, sparsity <= 1)
  switch(intdistr,
         normal = {
           stopifnot(length(intmean) == 1, length(intsd) == 1, intsd >= 0)
           intcoefs <- rnorm(n = nspecies^2, mean = intmean, sd = intsd)
         },
         uniform = {
           stopifnot(length(intrange) == 2, intrange[1] <= intrange[2])
           intcoefs <- runif(n = nspecies^2, min = intrange[1],
                             max = intrange[2])
         },
         {
           stop("'intdistr' should be 'normal' or 'uniform'.")
           intcoefs <- NULL
         }
  )
  intmat <- matrix(intcoefs, nrow = nspecies, ncol = nspecies)
  
  switch(selfintdistr,
         normal = {
           stopifnot(length(selfintmean) == 1, length(selfintsd) == 1,
                     selfintsd >= 0)
           # Draw variates from a truncated normal distribution to ensure
           # negative self-interactions without having to redraw for positive
           # deviates, using rtnorm() from the package TruncatedNormal.
           diag(intmat) <- rtnorm(n = nspecies, mu = selfintmean, sd = selfintsd,
                                  lb = -Inf, ub = 0, method = "invtransfo")
         },
         uniform = {
           stopifnot(length(selfintrange) == 2, selfintrange[1] <= selfintrange[2],
                     selfintrange[2] <= 0)
           diag(intmat) <- runif(n = nspecies,
                                 min = selfintrange[1], max = selfintrange[2])
         },
         {
           stop("'selfintdistr' should be 'normal' or 'uniform'.")
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
    sparse_index <- sample(1:(dim(indexmat)[1]), nsparseint)
    # Only keep rows of indexmat that were drawn from the sample. Use drop =
    # FALSE to prevent sparse_index becoming a vector of length two such that
    # two elements are set to zero if nsparseint == 1. 
    indexmat <- indexmat[sparse_index, , drop = FALSE]
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
#   - conjrate: numeric vector of length nspecies giving the intraspecies
#     conjugation rates for each species.
#   - taxmat: matrix where element in column n and row r gives the taxonomic
#     relatedness between donor species n and recipient species r. Although the
#     matrix is symmetric, this ordering is needed to correctly calculate the
#     conjugation rates.
#     Valid values are the character strings "SameSpecies", "SameFamily",
#     "SameOrder", "SameClass", and "OtherClass". Currently only "SameSpecies"
#     and, if taxmattypeset contains "PInOtherClass", "OtherClass" are used.
#     Interspecies conjugation rates are calculated based on the intraspecies
#     conjugation rates provided in conjrate and the taxonomic relatedness as
#     provided in this matrix taxmat.
#     By default R fills matrices by column, so byrow = TRUE should be used to
#     obtain the same matrix as the 'folded' vector.
#     An example matrix with E. coli, Klebsiella, and two Erwinia species:
#     taxmat <- matrix(c("SameSpecies", "SameFamily",  "SameOrder",   "SameOrder",
#                        "SameFamily",  "SameSpecies", "SameOrder",   "SameOrder",
#                        "SameOrder",   "SameOrder",   "SameSpecies", "SameSpecies",
#                        "SameOrder",   "SameOrder",   "SameSpecies", "SameSpecies"),
#                      nrow = 4, ncol = 4, byrow = TRUE)
#     Only the submatrix taxmat[1:nspecies, 1:nspecies] is used if nspecies <
#     max(nspeciesset). The use of multiple values for taxmat is NOT supported
#     (yet).
# Return:
#   - conjmat: matrix where the element in column n and row r gives the
#     conjugation rate from species n to species r.
# Notes:
# - To obtain a matrix where all values are identical (i.e., where conjugation
#   is not dependent on relatedness), fill taxmat with all values "SameSpecies"
#   and use all identical values for conjrate.
getconjmat <- function(nspecies, conjrate, taxmat) {
  stopifnot(all(diag(taxmat) == "SameSpecies"),
            isSymmetric.matrix(unname(taxmat)))
  conjratensp <- conjrate[1:nspecies]
  taxmatnsp <- taxmat[1:nspecies, 1:nspecies]
  
  # To obtain interspecies conjugation rates for the different levels of
  # taxonomic relatedness between donor and recipients, the intraspecies
  # conjugation rates are multiplied with the following conversion factors.
  rateconv <- c(SameSpecies = 1, SameFamily = 2.7, SameOrder = 0.1,
                SameClass = 0.05, OtherClass = 0.001)
  convmat <- matrix(NA, nrow = nspecies, ncol = nspecies)
  for(taxlevel in names(rateconv)) {
    convmat[which(taxmatnsp == taxlevel)] <- rateconv[taxlevel]
  }
  
  # Multiply column n of convmat giving the conversion factors for conjugation
  # from donor species n to the different recipient species, with element n of
  # conjrate. Using t(t(convmat) * conjrate) is faster for nspecies > 15
  conjmat <- convmat %*% diag(conjratensp)
  return(conjmat)
}

# Get equilibrium characteristics. On ecological and epidemiological stability
# see Roberts and Heesterbeek 2021.
geteqinfo <- function(model, abundance, intmat, growthrate,
                      cost = NULL, conjmat = NULL) {
  switch(model,
         "gLV" = {
           eigval <- eigen(x = jacobian.full(y = abundance, func = gLV,
                                             parms = list(growthrate = growthrate,
                                                          intmat = intmat)),
                           symmetric = FALSE, only.values = TRUE)$values
         },
         "gLVConj" = {
           eigval <- eigen(x = jacobian.full(y = abundance, func = gLVConj,
                                             parms = list(growthrate = growthrate,
                                                          intmat = intmat,
                                                          cost = cost,
                                                          conjmat = conjmat)),
                           symmetric = FALSE, only.values = TRUE)$values
         },
         "ecol" = {
           jac <- jacobian.full(y = abundance, func = gLVConj,
                                parms = list(growthrate = growthrate,
                                             intmat = intmat,
                                             cost = cost,
                                             conjmat = conjmat))
           jaclow <- jac[(nspecies + 1):(2*nspecies), 1:nspecies]
           if(!isTRUE(all.equal(range(jaclow), c(0, 0), check.attributes = FALSE))) {
             warning(paste("Jacobian matrix does not contain a block of zeros",
                     "in the lower-left corner,\nso currently used determination",
                     "of ecological stability is invalid"))
             print(jaclow)
           }
           eigval <- eigen(x = jac[1:nspecies, 1:nspecies],
                           symmetric = FALSE, only.values = TRUE)$values
         },
         "epi" = {
           jac <- jacobian.full(y = abundance, func = gLVConj,
                                parms = list(growthrate = growthrate, intmat = intmat,
                                             cost = cost, conjmat = conjmat))
           jaclow <- jac[(nspecies + 1):(2*nspecies), 1:nspecies]
           if(!isTRUE(all.equal(range(jaclow), c(0, 0), check.attributes = FALSE))) {
             warning(paste("Jacobian matrix does not contain a block of zeros",
                           "in the lower-left corner,\nso currently used",
                           "determination of epidemiological stability is invalid"))
             print(jaclow)
           }
           eigval <- eigen(x = jac[(nspecies + 1):(2*nspecies),
                                   (nspecies + 1):(2*nspecies)],
                           symmetric = FALSE, only.values = TRUE)$values
         },
         {
           eigval <- NULL
           warning("'model' should be one of 'gLV', 'gLVConj', 'ecol' or 'epi'")
         }
  )
  
  eigvalRep <- any(duplicated(eigval))
  # Using sort(eigval) to get eigenvalue with largest real part, because
  # max(eigval) does not work if eigval is complex. Complex values are sorted
  # first by the real part, then the imaginary part.
  eigval <- sort(eigval, decreasing = TRUE)[1]
  eigvalRe <- Re(eigval)
  eigvalIm <- Im(eigval)
  
  eqinfo <- c(eigvalRe = eigvalRe, eigvalIm = eigvalIm,
              eigvalReSign = sign(eigvalRe), eigvalImSign = sign(eigvalIm),
              eigvalRep = eigvalRep)
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
# pertpop is a character vector with the names of the populations to be
#   perturbed, e.g., pertpop = "R1" or pertpop = c("R1", "P1").
# pertmagn gives the absolute increase in populations for the perturbation
# tmax and tstep give the timesteps at which abundances should be calculated
#   (since variable step-size methods are used, those are not the only times
#   that integration occurs)
# If showplot == TRUE, the result is plotted (which is slow) or saved to a .png
#   file if saveplots == TRUE.
# If plotepistabwithP == TRUE, plots over time are plotted or saved to a
#   .png file if the equilibrium is epidemiologically stable but the total
#   number of plasmid-bearing bacteria at the end of the simulation is larger
#   than finalsmallstate.
# If verbose is TRUE, abundances before and after perturbation, and their
#   differences, are printed.
# Abundances grow to infinity when positive feedback is present, because the
#   system does not have an inherent carrying capacity. To prevent fatal errors
#   during the integration caused by state variables going to +Inf, or stepsizes
#   to 0, a rootfunction is used to terminate the integration if abundances get
#   very large. A warning is issued and the variable 'infgrowth' is set to 1 if
#   this occurs.
# When the simulation is finished, abundances smaller than finalsmallstate with
#   a negative derivative (i.e., small abundances that are declining) are set to
#   zero
# Output:
# Returns abunfinal, a list containing:
#   - R with absolute abundances of plasmid-free populations for each species,
#   - Rtotal with the sum of abundances of plasmid-free populations
#   - P with absolute abundances of plasmid-bearing population for each species
#     (or NULL if model == "gLV")
#   - Ptotal with the sum of abundances of plasmid-bearing populations (or NULL
#     if model == "gLV")
#   - timefinal indicating the last recorded time
#   - tmaxshort indicating tmax was not reached but no infinite growth occured
#   - eqreached indicating if equilibrium has been reached
#   - infgrowth indicating if infinite growth was detected
# If showplot == TRUE, plots over time are either shown or saved, depending on
#   the global variable saveplots
# Notes
#   Although this function returns absolute abundances, various of these are
#   converted to relative abundances later on in the script.
perturbequilibrium <- function(abundance, intmat, growthrate, cost, conjmat,
                               model, pertpop, pertmagn = 1000,
                               tmax = 1e4, tstep = 1, showplot = TRUE,
                               plotepistabwithP = FALSE, verbose = FALSE,
                               silentinfgrowth = FALSE,
                               silenteqnotreached = FALSE) {
  
  # Name abundances, set line type and colors, get derivatives of initial state
  # in the plasmid-free model.
  if(model == "gLV") {
    nspecies <- length(abundance)
    indexR <- 1:nspecies
    names(abundance) <- paste0(rep("R", nspecies), 1:nspecies)
    lty <- 1
    col <- mycol[1:nspecies]
    derivatives <- unlist(gLV(t = 0, n = abundance,
                              parms = list(growthrate = growthrate, intmat = intmat)))
  }
  
  if(model == "gLVConj") {
    nspecies <- length(abundance)/2
    indexR <- 1:nspecies
    indexP <- (nspecies + 1):(2*nspecies)
    names(abundance) <- c(paste0(rep("R", nspecies), 1:nspecies),
                          paste0(rep("P", nspecies), 1:nspecies))
    lty <- rep(c(1, 2), each = nspecies)
    col <- rep(mycol[1:nspecies], 2)
    if(!all(near(abundance[indexP], 0))) {
      warning("Initial state is NOT plasmid-free.")
    }
    # Derivatives of model without plasmids, to check if initial state is an
    # equilibrium in that model
    derivatives <- unlist(
      gLV(t = 0, n = abundance[indexR],
          parms = list(growthrate = growthrate, intmat = intmat))
    )
  }
  
  # Warn if initial plasmid-free state is not an equilibrium in the model
  # without plasmids
  if(!all(near(derivatives, 0))) {
    warntext <- paste("Initial (unperturbed) state is NOT an equilibrium in the
                      plasmid-free model! Derivatives are:",
                      paste(signif(derivatives), collapse = ", "))
    warning(warntext)
  }
  
  # Only perturb populations that exist
  if(!all(pertpop %in% names(abundance))) { 
    warning("Neglecting non-existent population(s) specified to perturb!")
    pertpop <- pertpop[which(pertpop %in% names(abundance))]
  }
  
  # Create perturbed abundances. Note that plasmid-baering bacteria replace
  # plasmid-free bacteria, instead of being added to the plasmid-free equilibrium.
  # In this way the number of bacteria for each species the same, such that the
  # ecological equilibrium is only affected by plasmid costs altering growth
  # rates, not by changes in abundances.
  abuninit <- abundance
  abunpert <- abundance
  abunpert[pertpop] <- abunpert[pertpop] + pertmagn
  pertpopP <- which(substr(pertpop, start = 1, stop = 1) == "P")
  # pertpopP becomes integer(0) in case of perturbation by plasmid-free bacteria
  if(length(pertpopP) > 0) {
    pertpopminus <- paste0("R", pertpopP)
    abunpert[pertpopminus] <- abunpert[pertpopminus] - pertmagn
  }
  
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
  tmaxshort <- 0
  if(!is.null(attributes(out)$troot)) {
    # One or more roots found
    
    if(any(attributes(out)$indroot == 1)) {
      eqreached <- 1  
    }
    
    if(any(attributes(out)$indroot == 2)) {
      infgrowth <- 1
      
      if(silentinfgrowth != TRUE) {
        warning(
          paste0("Integration was terminated at time = ",
                 round(attributes(out)$troot[which(attributes(out)$indroot == 2)], 2),
                 ", when the sum of abundances became",
                 "\nlarger than 1e10 times the initial total abundance,",
                 "indicating unbounded growth.\nAbundances then were ",
                 paste(names(abunfinaltemp), "=", signif(abunfinaltemp), collapse = ", "))
        )
      }
    }
  }
  
  if(eqreached == 0 && infgrowth != 1) {
    tmaxshort <- 1
    if(silenteqnotreached == FALSE) {
      warning("Equilibrium not reached. Increase tmax to prevent this? Final abundances were\n",
              paste(names(abunfinaltemp), "=", signif(abunfinaltemp), collapse = ", "))
    }
  }
  
  if(model == "gLV") {
    if(eqreached == 0) {
      # It does not make sense to store abundances in case of infinite
      # growth or if equilibrium is not reached, so record those as NA
      abunfinaltemp <- rep(NA, nspecies)
      names(abunfinaltemp) <- names(abundance)
      npopR <- NA
    } else {
      derivatives <- unlist(
        gLV(t = 0, n = abunfinaltemp,
            parms = list(growthrate = growthrate, intmat = intmat))
      )
      abunfinaltemp[which(derivatives < 0 & abunfinaltemp < finalsmallstate)] <- 0
      npopR <- length(which(abunfinaltemp > smallstate))
    }
    
    abunfinal <- list(R = abunfinaltemp, Rtotal = sum(abunfinaltemp),
                      npopR = npopR,
                      P = NULL, Ptotal = NULL, npopP = NULL,
                      timefinal = timefinal, tmaxshort = tmaxshort,
                      eqreached = eqreached, infgrowth = infgrowth)
  } else {
    if(eqreached == 0) {
      # It does not make sense to store abundances in case of infinite
      # growth or if equilibrium is not reached, so record those as NA
      abunfinaltemp <- rep(NA, 2*nspecies)
      names(abunfinaltemp) <- names(abundance)
      npopR <- NA
      npopP <- NA
    } else {
      derivatives <- unlist(
        gLVConj(t = 0, n = abunfinaltemp,
                parms = list(growthrate = growthrate, intmat = intmat,
                             cost = cost, conjmat = conjmat))
      )
      # Plasmid-bearing populations can reach densities smaller than smallstate
      # during the integration, because they are constantly produced through
      # conjugation after being set to 0.
      abunfinaltemp[which(derivatives < 0 & abunfinaltemp < finalsmallstate)] <- 0
      npopR <- length(which(abunfinaltemp[indexR] > smallstate))
      npopP <- length(which(abunfinaltemp[indexP] > smallstate))
    }
    
    abunfinal <- list(R = abunfinaltemp[indexR],
                      Rtotal = sum(abunfinaltemp[indexR]),
                      npopR = npopR,
                      P = abunfinaltemp[indexP],
                      Ptotal = sum(abunfinaltemp[indexP]),
                      npopP = npopP,
                      timefinal = timefinal, tmaxshort = tmaxshort,
                      eqreached = eqreached, infgrowth = infgrowth)
  }

  # eqinfoepi is created outside this function before this function is called.
  if(showplot == TRUE ||
     (plotepistabwithP == TRUE && model == "gLVConj" &&
      !is.na(abunfinal["Ptotal"]) && abunfinal["Ptotal"] > finalsmallstate &&
      eqinfoepi["eigvalRe"] < 0)) {
    subtitle <- paste0(abunmodel, ", intmean=", intmean, ",
                       selfintmean=", selfintmean, ", cost=", cost, ",
                       conjratecode=", conjratecode, ",\niter=", iter)
    if(saveplots == TRUE) {
      filename <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M_%OS3"), ".png")
      png(filename = filename, height = 785)
    }
    matplot.deSolve(out, lty = lty, col = col, ylab = "Abundance",
                    log = "y", sub = subtitle, lwd = 2, legend = list(x = "bottomright"))
    grid()
    abline(h = abuninit)
    if(saveplots == TRUE) {
      dev.off()
    }
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
  
  return(abunfinal)
}

# List of functions called within dplyr::summarise() to get summary statistics
# for the input data. If all values in x are NA (which is the case for example
# when data on species 3 is considered in the two-species model), NA will be
# returned, otherwise the NA values are dropped and the function will be applied
# to the remaining values.
getsummary4 <- list(
  min = ~if(any(!is.na(.x))) {min(.x, na.rm = TRUE)} else {NA},
  mean = ~if(any(!is.na(.x))) {mean(.x, na.rm = TRUE)} else {NA},
  median = ~if(any(!is.na(.x))) {median(.x, na.rm = TRUE)} else {NA},
  max = ~if(any(!is.na(.x))) {max(.x, na.rm = TRUE)} else {NA}
)
getfracnotzero <- list(frac = ~mean(.x != 0))

# Function to create plots. The plotted object is returned, such that it can be
# further modified like any other ggplot object.
CreatePlot <- function(dataplot = plotdata, xvar = "intmean", yvar = "selfintmean",
                       contour_var = NULL, contour_breaks = 0.5, contour_col = NULL,
                       contour_lty = NULL,
                       limits = NULL, limx = NULL, limy = NULL, ratio = 1,
                       fillvar, filltitle, filltype = "discrete",
                       title = NULL, subtitle = NULL,
                       labx = "Mean interaction coefficient",
                       laby = "Mean selfinteraction coefficient",
                       tag = NULL, addstamp = FALSE, diagonal = "none",
                       linezero = TRUE,
                       facetx = "taxmatcode + abunmodelcode + cost",
                       facety = "conjratecode + nspecies",
                       dropfacets = TRUE,
                       as.table = TRUE,
                       marginx = NULL, marginy = NULL, base_size = 11,
                       rotate_x_labels = TRUE, rotate_legend = FALSE,
                       save = saveplots, filename = NULL, ...) {
  caption <- paste(unique(dataplot$niter), "iterations")
  if(exists("DateTimeStamp") == FALSE) {
    DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
    if(addstamp == TRUE) {
      warning(paste("DateTimeStamp created to include in plot does not",
                    "correspond to filename of the dataset"))
    }
  }
  if(addstamp == TRUE) {
    caption <- paste0(caption, ", ", DateTimeStamp)
  }
  
  # Facets that have only a single unique value are not shown if dropfacets = TRUE
  if(dropfacets == TRUE) {
    if(!is.null(facetx) && facetx != ".") {
      elements_facetx <- unlist(strsplit(facetx, split = " + ", fixed = TRUE))
      if(!all(elements_facetx %in% colnames(dataplot))) {
        warning("Dropped facetx variables that are not column names of plotdata")
        elements_facetx <- elements_facetx[elements_facetx %in% colnames(dataplot)]
      }
      include_facetx <- rep(FALSE, length(elements_facetx))
      for(index_facetx in seq_along(elements_facetx)) {
        if(length(unique(dataplot[, elements_facetx[index_facetx]])) > 1) {
          include_facetx[index_facetx] <- TRUE
        }
      }
      facetx <- paste0(elements_facetx[include_facetx], collapse = " + ")
      if(facetx == "") {
        facetx <- "."
      }
    }
    
    if(!is.null(facety) && facety != ".") {
      elements_facety <- unlist(strsplit(facety, split = " + ", fixed = TRUE))
      if(!all(elements_facety %in% colnames(dataplot))) {
        warning("Dropped facety variables that are not column names of plotdata")
        elements_facety <- elements_facety[elements_facety %in% colnames(dataplot)]
      }
      include_facety <- rep(FALSE, length(elements_facety))
      for(index_facety in seq_along(elements_facety)) {
        if(length(unique(dataplot[, elements_facety[index_facety]])) > 1) {
          include_facety[index_facety] <- TRUE
        }
      }
      facety <- paste0(elements_facety[include_facety], collapse = " + ")
      if(facety == "") {
        facety <- "."
      }
    }
  }
  
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar)) + 
    theme_bw(base_size = base_size) +
    facet_grid(as.formula(paste(facety, "~", facetx)), as.table = as.table,
               labeller = mylabeller) +
    theme(legend.position = "bottom") +
    labs(title = title, subtitle = subtitle,
         x = labx, y = laby, caption = caption, tag = tag)
  if(!is.null(limx)) {
    p <- p + scale_x_continuous(limits = limx, expand = c(0, 0))
  } else {
    p <- p + scale_x_continuous(expand = c(0, 0))
  }
  if(!is.null(limy)) {
    p <- p + scale_y_continuous(limits = limy, expand = c(0, 0))
  } else {
    p <- p + scale_y_continuous(expand = c(0, 0))
  }
  if(!is.null(ratio)) {
    p <- p + coord_fixed(ratio = ratio, expand = FALSE)
  }
  if(!is.null(marginx)) {
    p <- p + theme(strip.text.x = element_text(margin = margin(marginx)))
  }
  if(!is.null(marginy)) {
    p <- p + theme(strip.text.y = element_text(margin = margin(marginy)))
  }
  if(!is.null(contour_var)) {
    p <- p + geom_contour(aes_string(z = contour_var, color = contour_col,
                                     linetype = contour_lty),
                          breaks = contour_breaks, size = 1)
  } else {
    p <- p + geom_raster()
    if(filltype == "discrete") {
      p <- p + scale_fill_viridis_d(filltitle,
                                    limits = if(is.null(limits)) {
                                      as.factor(c(-1, 1))
                                    } else {
                                      factor(limits)
                                    },
                                    labels = filllabels)
    }
    if(filltype == "binned") {
      p <- p + scale_fill_viridis_b(filltitle, breaks = limits)
    }
    if(filltype == "continuous") {
      p <- p + scale_fill_viridis_c(filltitle, limits = limits)
    }
  }
  if(rotate_x_labels == TRUE) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
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
  if(linezero == TRUE) {
    p <- p + geom_vline(xintercept = 0, col = "white", size = 1.1)
  }
  if(save == TRUE) {
    if(is.null(filename)) {
      filename <- gsub("/", ".", fillvar)
      filename <- gsub("\\*", ".", filename)
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
  return(p)
}


#### Testing functions ####
# 
# nspecies <- 4
# abunbrokenstick <- brokenstick(nspecies = nspecies, totalabun = totalabun,
#                                takelimit = TRUE)
# abundompreempt <- dompreempt(nspecies = nspecies, totalabun = totalabun,
#                              takelimit = TRUE)
# abunbrokenstick
# abundompreempt
# 
# (intmat1 <- getintmat(nspecies = nspecies))
# (growthratebrokenstick <- getgrowthrate(abundance = abunbrokenstick, intmat = intmat1))
# (growthratedompreempt <- getgrowthrate(abundance = abundompreempt, intmat = intmat1))
# checkequilibrium(abundance = abunbrokenstick, intmat = intmat1,
#                  growthrate = growthratebrokenstick, printderivatives = TRUE,
#                  showplot = TRUE) # At equilibrium
# checkequilibrium(abundance = 0.9*abunbrokenstick, intmat = intmat1,
#                  growthrate = growthratebrokenstick, printderivatives = TRUE,
#                  showplot = TRUE) # Not at equilibrium

# geteqinfo returns eigvalRe, eigvalIm, eigvalReSign, eigvalImSign, and eigvalRep
# geteqinfo(model = "gLV", abundance = abunbrokenstick, intmat = intmat1,
#           growthrate = growthratebrokenstick)
# geteqinfo(model = "gLV", abundance = abundompreempt, intmat = intmat1,
#           growthrate = growthratedompreempt)


#### Running the simulations ####
set.seed(seed = 314, kind = "default", normal.kind = "default", sample.kind = "default")
starttime <- Sys.time()

# Create matrix to store data
nrowplotdata <- prod(lengths(list(nspeciesset, abunmodelset, intmeanset,
                                  selfintmeanset, costset, conjrateset, taxmattypeset),
                             use.names = FALSE))
ncolplotdata <- if(simulateinvasion == TRUE) {
  17*4 + 4*4*maxnspecies + 37
} else {
  3*4 + 29
}
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = ncolplotdata)
nrowdatatotal <- prod(lengths(list(abunmodelset,intmeanset, selfintmeanset,
                                   costset, conjrateset, taxmattypeset),
                              use.names = FALSE))*niter*sum(nspeciesset)
datatotal <- matrix(data = NA, nrow = nrowdatatotal, ncol = 57 + 4*maxnspecies)
indexdatatotal <- 1

# Run simulations
rowindexplotdata <- 1
rowindexdata <- 1
for(nspecies in nspeciesset) {
  if(PInMostAbun == TRUE) {
    pertpop <- "R1"
    pertpopconj <- "P1"
  } else {
    pertpop <- paste0("R", nspecies)
    pertpopconj <- paste0("P", nspecies)
  }

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
      # print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
      #              ", intmean = ", intmean,
      #              ": started at ", Sys.time()), quote = FALSE)
      
      for(selfintmean in selfintmeanset) {
        print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
                     ", intmean = ", intmean, ", selfintmean = ", selfintmean,
                     ": started at ", Sys.time()), quote = FALSE)
        nrowdata <- niter * nspecies * length(costset) * length(conjrateset) *
          length(taxmattypeset)
        data <- matrix(data = NA, nrow = nrowdata, ncol = 57 + 4*maxnspecies)
        indexdata <- 1
        relabunRsp <- rep(NA, maxnspecies)
        relabunRconjsp <- rep(NA, maxnspecies)
        relabunPconjsp <- rep(NA, maxnspecies)

        for(iter in 1:niter) {
          stableeq <- FALSE
          iterintmat <- 0
          conjratecode <- NA
          
          # If niterintmat > 1 (which is not the default), niterintmat attempts
          # are done to get a combination of interaction matrix and growth rates
          # that results in a stable plasmid-free equilibrium in the model
          # without conjugation.
          while(stableeq == FALSE && iterintmat < niterintmat) {
            intmat <- getintmat(nspecies = nspecies,
                                intmean = intmean, selfintmean = selfintmean)
            growthrate <- getgrowthrate(abundance = abundance, intmat = intmat)
            eqinfo <- geteqinfo(model = "gLV", abundance = abundance,
                                intmat = intmat, growthrate = growthrate)
            if(eqinfo["eigvalRe"] < 0) {
              stableeq <- TRUE
            }
            iterintmat <- iterintmat + 1
          }
          # To get a warning if the plasmid-free equilibrium is not stable,
          # uncomment next lines.
          # if(stableeq == FALSE) {
          # warning(paste("No stable equilibrium has been found in",
          #               niterintmat, "attempts."))
          # }
          
          # Simulate invasion with plasmid-free bacteria of plasmid-free
          # equilibrium in model without conjugation (i.e., test internal
          # stability)
          if(simulateinvasion == TRUE) {
            if(stableeq == FALSE) {
              abunfinal <- perturbequilibrium(abundance = abundance, intmat = intmat,
                                              growthrate = growthrate,
                                              cost = NULL, conjmat = NULL,
                                              model = "gLV", pertpop = pertpop,
                                              pertmagn = 1, tmax = 1e4,
                                              showplot = FALSE, verbose = FALSE,
                                              silentinfgrowth = TRUE,
                                              silenteqnotreached = TRUE)
            } else {
              # No need for simulations if equilibrium is stable
              abunfinal <- list(R = abundance, Rtotal = sum(abundance),
                                npopR = nspecies,
                                P = NULL, Ptotal = NULL, npopP = NULL,
                                timefinal = 1, tmaxshort = 0,
                                eqreached = 1, infgrowth = 0)
            }
          } else {
            # No simulations over time performed, so set values to NA
            abunfinal <- list(R = rep(NA, nspecies), Rtotal = NA,
                              npopR = NA,
                              P = NULL, Ptotal = NULL, npopP = NULL,
                              timefinal = NA, tmaxshort = NA,
                              eqreached = NA, infgrowth = NA)
          }
          
          # Model without conjugation, so Rtotal is total abundance as P does
          # not exist 
          relabunRsp[1:nspecies] <- abunfinal$R / abunfinal$Rtotal

          for(cost in costset) {
            conjratecode <- 0
            
            for(conjrate in conjrateset) {
              conjratecode <- conjratecode + 1
              taxmatcode <- 0
              
            for(taxmattype in taxmattypeset) {
                taxmatcode <- taxmatcode + 1
                taxmat <- matrix(rep("SameSpecies", nspecies^2),
                                 nrow = nspecies, ncol = nspecies, byrow = TRUE)
                if(taxmattype == "PInOtherClass") {
                  if(PInMostAbun == TRUE) {
                    # The plasmid is introduced in the most abundant species,
                    # so set the first row and column of the matrix to OtherClass,
                    # with exception of the diagonal which should remain SameSpecies
                    taxmat[1, -1] <- "OtherClass"
                    taxmat[-1, 1] <- "OtherClass"                    
                  } else {
                    # The plasmid is introduced in the least abundant species,
                    # so set the last row and column of the matrix to OtherClass,
                    # with exception of the diagonal which should remain SameSpecies
                    taxmat[nspecies, -nspecies] <- "OtherClass"
                    taxmat[-nspecies, nspecies] <- "OtherClass"
                  }
                }
              conjmat <- getconjmat(nspecies = nspecies,
                                    conjrate = conjrate, taxmat = taxmat)
              
              # Get equilibrium characteristics for plasmid-free equilibrium in
              # the model with conjugation
              eqinfoconj <- geteqinfo(model = "gLVConj",
                                      abundance = c(abundance, rep(0, nspecies)),
                                      intmat = intmat, growthrate = growthrate,
                                      cost = cost, conjmat = conjmat)
              
              # Get equilibrium characteristics regarding ecological and
              # epidemiological stability of the plasmid-free equilibrium (see
              # Roberts and Heesterbeek 2021).
              eqinfoecol <- geteqinfo(model = "ecol",
                                      abundance = c(abundance, rep(0, nspecies)),
                                      intmat = intmat, growthrate = growthrate,
                                      cost = cost, conjmat = conjmat)
              eqinfoepi <- geteqinfo(model = "epi",
                                     abundance = c(abundance, rep(0, nspecies)),
                                     intmat = intmat, growthrate = growthrate,
                                     cost = cost, conjmat = conjmat)

              # Append the signs of the real parts of the largest eigenvalues to
              # indicate change in stability without and with the plasmid (-1 =
              # stable, 0 = neutral, 1 = unstable). The + 3 ensure numbers are
              # positive, to prevent paste0(-1, -1) leading to NA because of
              # as.integer(-1-1))
              compstability <- as.integer(paste0(eqinfo["eigvalReSign"] + 3,
                                                 eqinfoconj["eigvalReSign"] + 3))
              
              # To simulate invasion of the plasmid-free equilibrium with plasmids,
              # the abundances of the plasmid-free populations have to be appended
              # to the abundances of the plasmid-bearing populations
              if(simulateinvasion == TRUE) {
                if(eqinfoconj["eigvalRe"] >= 0) {
                  abunfinalconj <- perturbequilibrium(abundance = c(abundance, rep(0, nspecies)),
                                                      intmat = intmat, growthrate = growthrate,
                                                      cost = cost, conjmat = conjmat,
                                                      model = "gLVConj", pertpop = pertpopconj,
                                                      pertmagn = 1, tmax = 1e4, 
                                                      showplot = FALSE, verbose = FALSE,
                                                      silentinfgrowth = TRUE,
                                                      silenteqnotreached = TRUE)
                } else {
                  # No need for simulations if equilibrium is stable
                  abunfinalconj <- list(R = abundance, Rtotal = sum(abundance),
                                        npopR = nspecies,
                                        P = rep(0, nspecies), Ptotal = 0,
                                        npopP = 0,
                                        timefinal = 1, tmaxshort = 0,
                                        eqreached = 1, infgrowth = 0)
                }
                
              } else {
                # No simulations over time performed, so set values to NA
                abunfinalconj <- list(R = rep(NA, nspecies), Rtotal = NA,
                                  npopR = NA,
                                  P = rep(NA, nspecies), Ptotal = NA,
                                  npopP = NA,
                                  timefinal = NA, tmaxshort = NA,
                                  eqreached = NA, infgrowth = NA)
              }
              
              # Total (cells / mL) and relative (fractions of total) abundances
              abuntotalconj <- abunfinalconj$Rtotal + abunfinalconj$Ptotal
              relabunRconjsp[1:nspecies] <- abunfinalconj$R / abuntotalconj
              relabunPconjsp[1:nspecies] <- abunfinalconj$P / abuntotalconj

              # Using abunfinalconj$R + abunfinalconj$P > smallstate if the
              # abundances are NA leads to 0 instead of NA, so instead use
              # if-else construct.
              if(abunfinalconj$eqreached == 0 | simulateinvasion == FALSE) {
                nspeciesconj <- NA 
                fracPformedbypertpop <- NA
              } else {
                nspeciesconj <- length(which(
                  abunfinalconj$R + abunfinalconj$P > smallstate
                ))
                if(abunfinalconj$Ptotal > 0) {
                  fracPformedbypertpop <- sum(abunfinalconj$P[pertpopconj]) /
                    abunfinalconj$Ptotal
                } else {
                  fracPformedbypertpop <- NA
                }
              }
              
              indexdatanew <- indexdata + nspecies
              
              data[indexdata:(indexdatanew - 1), ] <- cbind(
                niter, nspecies, abunmodelcode, intmean, selfintmean,
                cost, conjratecode, taxmatcode, iter, 1:nspecies, abundance,
                diag(intmat), c(growthrate), iterintmat,
                matrix(rep(eqinfo, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(eqinfoconj, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(eqinfoecol, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(eqinfoepi, nspecies), nrow = nspecies, byrow = TRUE),
                compstability,
                abunfinal$infgrowth, abunfinalconj$infgrowth,
                abunfinal$eqreached, abunfinalconj$eqreached,
                abunfinal$tmaxshort, abunfinalconj$tmaxshort,
                abunfinal$timefinal, abunfinalconj$timefinal,
                abunfinal$Rtotal, abunfinalconj$Rtotal, abunfinalconj$Ptotal,
                abuntotalconj, abunfinalconj$Rtotal/abuntotalconj,
                abunfinal$npopR, abunfinalconj$npopR, abunfinalconj$npopP,
                abunfinalconj$npopR + abunfinalconj$npopP,
                abunfinalconj$Ptotal / abuntotalconj,
                nspeciesconj, abunfinalconj$npopR / nspeciesconj,
                abunfinalconj$npopP / nspeciesconj,
                fracPformedbypertpop,
                matrix(rep(relabunRsp, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(relabunRconjsp, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(relabunPconjsp, nspecies), nrow = nspecies, byrow = TRUE),
                matrix(rep(relabunRconjsp + relabunPconjsp, nspecies),
                       nrow = nspecies, byrow = TRUE)
              )
              indexdata <- indexdatanew
            }
            }
          }
        }
        datatotal[indexdatatotal:(indexdatatotal + nrowdata - 1), ] <- data
        indexdatatotal <- indexdatatotal + nrowdata
        
        colnames(data) <- c("niter", "nspecies", "abunmodelcode",
                            "intmean", "selfintmean", "cost", "conjratecode",
                            "taxmatcode", "iter", "species", "abundance",
                            "selfintdata", "growthrate",
                            "iterintmat",
                            "eigvalRe", "eigvalIm", "eigvalReSign",
                            "eigvalImSign", "eigvalRep",
                            "eigvalReconj", "eigvalImconj",
                            "eigvalReSignconj", "eigvalImSignconj", "eigvalRepconj",
                            "eigvalReecol", "eigvalImecol",
                            "eigvalReSignecol", "eigvalImSignecol", "eigvalRepecol",
                            "eigvalReepi", "eigvalImepi",
                            "eigvalReSignepi", "eigvalImSignepi", "eigvalRepepi",
                            "compstability",
                            "infgrowth", "infgrowthconj",
                            "eqreached", "eqreachedconj",
                            "tmaxshort", "tmaxshortconj",
                            "timefinal", "timefinalconj",
                            "abuntotal", "abunRtotalconj", "abunPtotalconj",
                            "abuntotalconj", "relabunRconj",
                            "npopR", "npopRconj", "npopPconj", "npopconj",
                            "relabunPconj", "nspeciesconj",
                            "fracspeciesRconj", "fracspeciesPconj",
                            "fracPformedbypertpop",
                            paste0("relabunRsp", 1:maxnspecies),
                            paste0("relabunRconjsp", 1:maxnspecies),
                            paste0("relabunPconjsp", 1:maxnspecies),
                            paste0("relabunconjsp", 1:maxnspecies))

        # Get summary data which do not depend on simulated invasion for all
        # combinations of costs and conjugation rates
        summarydata <- as_tibble(data) %>%
          group_by(cost, conjratecode, taxmatcode) %>%
          summarise(
            across(c(selfintdata, growthrate, iterintmat),
                   getsummary4, .names = "{.col}{.fn}"),
            fracstable = mean(eigvalRe < 0),
            fracstableconj = mean(eigvalReconj < 0),
            fracstableecol = mean(eigvalReecol < 0),
            fracstableepi = mean(eigvalReepi < 0),
            fracreal = mean(eigvalIm == 0),
            fracrealconj = mean(eigvalImconj == 0),
            fracrealecol = mean(eigvalImecol == 0),
            fracrealepi = mean(eigvalImepi == 0),
            across(c(eigvalRep, eigvalRepconj, eigvalRepecol, eigvalRepepi),
                   getfracnotzero, .names = "{.fn}{.col}"),
            fracstablestable = mean(compstability == 22),
            fracstableneutral = mean(compstability == 23),
            fracstableunstable = mean(compstability == 24),
            fracneutralstable = mean(compstability == 32),
            fracneutralneutral = mean(compstability == 33),
            fracneutralunstable = mean(compstability == 34),
            fracunstablestable = mean(compstability == 42),
            fracunstableneutral = mean(compstability == 43),
            fracunstableunstable = mean(compstability == 44),
            .groups = "drop"
          )
        
        if(simulateinvasion == TRUE) {
          # If invasion was simulated, get summary of generated data for all
          # combinations of costs and conjugation rates
          summarydatasimulation <- as_tibble(data) %>%
            group_by(cost, conjratecode, taxmatcode) %>%
            summarise(
              across(c(infgrowth, infgrowthconj, eqreached, eqreachedconj,
                       tmaxshort, tmaxshortconj),
                     getfracnotzero, .names = "{.col}{.fn}"),
              timefinalmedian = median(timefinal),
              timefinalconjmedian = median(timefinalconj),
              across(c(abuntotal, abunRtotalconj,
                       starts_with("relabunR"), abunPtotalconj,
                       starts_with("relabunP"), abuntotalconj,
                       relabunRconj, starts_with("relabunconjsp"),
                       starts_with("npop"), nspeciesconj,
                       fracspeciesRconj, fracspeciesPconj, fracPformedbypertpop),
                     getsummary4, .names = "{.col}{.fn}"),
              .groups = "drop"
            )
          summarydata <- full_join(summarydata, summarydatasimulation,
                                   by = c("cost", "conjratecode", "taxmatcode"))
        }
        
        rowindexplotdatanew <- rowindexplotdata + length(costset) *
          length(conjrateset) * length(taxmattypeset)
        plotdata[rowindexplotdata:(rowindexplotdatanew - 1), ] <- as.matrix.data.frame(
          tibble(niter, nspecies, abunmodelcode, intmean, selfintmean, summarydata))
        rowindexplotdata <- rowindexplotdatanew
      }
    }
  }
}
duration <- Sys.time() - starttime
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)

colnames(plotdata) <- c("niter", "nspecies", "abunmodelcode",
                        "intmean", "selfintmean", colnames(summarydata))
colnames(datatotal) <- colnames(data)
if(simulateinvasion == TRUE) {
  eqnotreached <- 1 - plotdata[, "eqreachedfrac"]
  eqnotreachedconj <- 1 - plotdata[, "eqreachedconjfrac"]
  if(any(eqnotreached > 0)) {
    warning(paste0("Fraction of iterations where equilibrium has not been reached ",
                   "after pertubation with plasmid-free\nbacteria ranges from ",
                   summary(eqnotreached)["Min."], " to ",
                   summary(eqnotreached)["Max."]," (mean: ",
                   summary(eqnotreached)["Mean"], "; median: ",
                   summary(eqnotreached)["Median"], "). ",
                   "Use silenteqnotreached = FALSE in perturbequilibrium() for more info"))
  }
  if(any(eqnotreachedconj > 0)) {
    warning(paste0("Fraction of iterations where equilibrium has not been reached ",
                   "after pertubation with plasmid-bearing\nbacteria ranges from ",
                   summary(eqnotreachedconj)["Min."], " to ",
                   summary(eqnotreachedconj)["Max."], " (mean: ",
                   summary(eqnotreachedconj)["Mean"], "; median: ",
                   summary(eqnotreachedconj)["Median"], "). ",
                   "Use silenteqnotreached = FALSE in perturbequilibrium() for more info"))
  }
}
warnings()


#### Saving output to .csv files ####
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
if(PInMostAbun == FALSE) {
  DateTimeStamp <- paste0(DateTimeStamp, "PInLeastAbun")
}

write.csv(plotdata, file = paste0(DateTimeStamp, "multispecies.csv"),
          quote = FALSE, row.names = FALSE)
if(nrow(datatotal) < 250000) {
  write.csv(datatotal, file = paste0(DateTimeStamp, "multispeciestotal.csv"),
            quote = FALSE, row.names = FALSE)
}

# Saving settings
names(conjrateset) <- paste0("conjrateset", 1:length(conjrateset))
settings <- c(list(niter = niter, niterintmat = niterintmat,
                   simulateinvasion = simulateinvasion,
                   smallstate = smallstate, smallchange = smallchange,
                   tstep = formals(perturbequilibrium)$tstep,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   costset = costset, conjrateset, taxmattype = taxmattypeset,
                   costtype = costtype, PInMostAbun = PInMostAbun,
                   duration = duration))
for(index in 1:length(settings)) {
  write.table(t(as.data.frame(settings[index])), 
              paste0(DateTimeStamp, "settings.csv"), append = TRUE,
              quote = FALSE, sep = ",", col.names = FALSE)
}


#### Reading previously saved data from a .csv-file ####
## To read data from csv-file, uncomment this section and fill in the 
## needed datetimestamp
# filename <- "2021_05_04_17_44multispecies.csv"
# plotdata <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                      dec = ".", stringsAsFactors = FALSE)
# # If plotdata has only one column, probably a semicolon instead of a comma is
# # used as separator in .csv-files. So read the file again.
# if(dim(plotdata)[2] == 1) {
#   plotdata <- read.csv(filename, header = TRUE, sep = ";", quote = "\"",
#                        dec = ".", stringsAsFactors = FALSE)
# }
# plotdata <- as.data.frame(plotdata)
# DateTimeStamp <- substr(filename, 1, 16)
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))


#### Labels and limits for plots ####
labspecies <- paste("Species", 1:maxnspecies)
names(labspecies) <- 1:maxnspecies
labnspecies <- paste(nspeciesset, "species")
names(labnspecies) <- nspeciesset
labmodel <- c("Broken stick", "Dom. preemption")
names(labmodel) <- c(1, 2)
labcost <- paste0("Cost: ", costset, "/h")
names(labcost) <- costset
labconjrate <- paste("Conjset", 1:length(conjrateset))
names(labconjrate) <- 1:length(conjrateset)
labtaxmat <- taxmattypeset
names(labtaxmat) <- 1:length(taxmattypeset)
mylabeller <- labeller(species = labspecies, nspecies = labnspecies,
                       abunmodelcode = labmodel,
                       cost = labcost, conjratecode = labconjrate,
                       taxmatcode = labtaxmat, .default = label_value)

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
# The error '$ operator is invalid for atomic vectors' arises if the matrix
# 'plotdata' has not been comverted to a dataframe. To convert it to a dataframe,
# run plotdata <- as.data.frame(plotdata)


# ## Show border of ecological stability with heatmap in CreatePlot()
# CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableecol",
#            filltitle = "fracstableecol", filltype = "continuous",
#            limy = c(1, length(conjrateset)), ratio = NULL,
#            labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
#            facetx = "taxmatcode + intmean", facety = "nspecies",
#            rotate_x_labels = TRUE,
#            filename = "ecostabxcostyconj")
# 
# ## Show border of ecological stability with a contour plot in CreatePlot()
# # Values in a facet are either all stable, or all unstable (see heatmap above),
# # making it impossible to plot contours delimiting stable and unstable regions.
# CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableecol",
#            contour_var = "fracstableecol", contour_col = "as.factor(nspecies)",
#            contour_lty = "as.factor(intmean)",
#            limy = c(1, length(conjrateset)), ratio = NULL,
#            labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
#            facetx = "taxmatcode", facety = "nspecies",
#            rotate_x_labels = FALSE, save = FALSE) +
#   guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1))


if(bifurparms == TRUE) {
  ## Show border of epidemiological stability with a contour plot in CreatePlot()
  # I set save to FALSE and used ggsave() to ensure the added guides arguments are
  # included in the saved plots.
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = "taxmatcode + intmean", facety = "nspecies",
             rotate_x_labels = TRUE, save = FALSE) +
    guides(col = guide_legend(ncol = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatintmeanynspecies.png"))
  }
  
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(intmean)",
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = "taxmatcode", facety = "nspecies",
             rotate_x_labels = FALSE, save = FALSE) +
    guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspecies.png"))
  }
  # So intmean does not affect the border of stability in conjugation rate/cost-space
  
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = "taxmatcode", facety = "intmean",
             rotate_x_labels = TRUE, save = FALSE) +
    guides(col = guide_legend(ncol = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatyintmean.png"))
  }
  
  # Invasion is possible for very slightly higher costs for a given conjugation
  # rate in a microbiome of 2 species than in a microbiome of 4 or 6 species when
  # the initially plasmid-bearing species belongs to another class.
  
  # When the plasmid-bearing species belongs to another class, invasion is only
  # possible if costs are slightly lower for a given conjugation rate than when
  # all populations belonging to the same species.
  
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(intmean)",
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = "taxmatcode", facety = ".",
             rotate_x_labels = FALSE, save = FALSE) +
    guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmat.png"))
  }
  
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(taxmatcode)",
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = ".", facety = "nspecies",
             rotate_x_labels = FALSE, save = FALSE) +
    guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabynspecies.png"))
  }
  
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(taxmatcode)",
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = ".", facety = ".",
             rotate_x_labels = FALSE, save = FALSE) +
    guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistab.png"))
  }
  
  # Need to set filltype to continuous to prevent error on missing filllabels
  CreatePlot(xvar = "cost", yvar = "conjratecode", fillvar = "fracstableepi",
             filltitle = "fracstableepi", contour_var = NULL, contour_col = NULL,
             contour_lty = NULL, filltype = "continuous", limx = range(costset),
             limy = c(1, length(conjrateset)), ratio = NULL,
             title = "Epidemiological stability",
             labx = "Cost", laby = "Conjugation rate code", linezero = FALSE,
             facetx = "taxmatcode", facety = "nspecies",
             rotate_x_labels = FALSE, filename = "epistabheatmap")
}


## Compare equilibrium characteristics for the models without and with plasmids.
# If invasion has been simulated, data on infinite growth and reaching
# equilibrium after perturbation is also plotted.
CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction,
           facetx = "abunmodelcode", facety = "nspecies",
           filename = "fracstablefewfacets")
CreatePlot(fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracstableconj",
           filltitle = "Fraction stable\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "1 - fracstable", filltitle = "Fraction unstable",
           filltype = "continuous", limits = limitsfraction,
           facetx = "abunmodelcode", facety = "nspecies",
           filename = "fracunstablefewfacets")
CreatePlot(fillvar = "1 - fracstable", filltitle = "Fraction unstable",
           filltype = "continuous", limits = limitsfraction,
           filename = "fracunstablecontinuous")
CreatePlot(fillvar = "1 - fracstableconj",
           filltitle = "Fraction unstable\nwith conjugation",
           filltype = "continuous", limits = limitsfraction,
           filename = "fracunstableconjcontinuous")

CreatePlot(fillvar = "fracstableecol",
           filltitle = "Fraction ecologically\nstable",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracstableepi",
           filltitle = "Fraction epidemiologically\nstable",
           filltype = "continuous", limits = limitsfraction)

# Show the effect of adding conjugation on stability
CreatePlot(fillvar = "fracunstableunstable + fracneutralneutral + fracstablestable",
           filltitle = "Fraction stability the same\nwithout and with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracunstableunstable",
           filltitle = "Fraction unstable without and\nunstable with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracunstableneutral",
           filltitle = "Fraction unstable without but\nneutral with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracunstablestable",
           filltitle = "Fraction unstable without but\nstable with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracneutralunstable",
           filltitle = "Fraction neutral without but\nunstable with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracneutralneutral",
           filltitle = "Fraction neutral without and\nneutral with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracneutralstable",
           filltitle = "Fraction neutral without but\nstable with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracstableunstable",
           filltitle = "Fraction stable without but\nunstable with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracstableneutral",
           filltitle = "Fraction stable without but\nneutral with conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracstablestable",
           filltitle = "Fraction stable without and\nstable with conjugation",
           filltype = "continuous", limits = limitsfraction)


# Show dominant eigenvalues
datatotalfilteredspecies <- filter(datatotal, near(species, 1))
datatotalfilteredspeciesiter <- filter(datatotalfilteredspecies, near(iter, 1))
limitseigvalRe <- range(c(datatotalfilteredspeciesiter[, "eigvalRe"],
                          datatotalfilteredspeciesiter[, "eigvalReconj"]))
limitseigvalIm <- range(c(datatotalfilteredspeciesiter[, "eigvalIm"],
                          datatotalfilteredspeciesiter[, "eigvalImconj"]))
CreatePlot(dataplot = datatotalfilteredspeciesiter,
           fillvar = "eigvalRe",
           filltitle = "Real part of\ndominant eigenvalue",
           filltype = "continuous", limits = limitseigvalRe)
CreatePlot(dataplot = datatotalfilteredspeciesiter,
           fillvar = "eigvalReconj",
           filltitle = "Real part of\ndominant eigenvalue\nwith conjugation",
           filltype = "continuous", limits = limitseigvalRe)

limitseigvalRe <- range(c(datatotalfilteredspecies[, "eigvalRe"],
                          datatotalfilteredspecies[, "eigvalReconj"]))
limitseigvalIm <- range(c(datatotalfilteredspecies[, "eigvalIm"],
                          datatotalfilteredspecies[, "eigvalImconj"]))
ggplot(data = datatotalfilteredspecies,
       aes(x = eigvalRe, y = eigvalIm,
           color = as.factor(eigvalReSign))) +
  scale_x_continuous(limits = limitseigvalRe) +
  scale_y_continuous(limits = limitseigvalIm) +
  geom_point() +
  theme(legend.position = "bottom") +
  labs(x = "Real part dominant eigenvalue",
       y = "Imaginary part dominant eigenvalue",
       caption = paste(niter, "iterations")) +
  facet_grid(rows = vars(nspecies, conjratecode),
             cols = vars(taxmatcode, abunmodelcode, cost),
             labeller = mylabeller)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "eigenvaluesdistr.png"))
}

ggplot(data = datatotalfilteredspecies,
       aes(x = eigvalReconj, y = eigvalImconj,
           color = as.factor(eigvalReSignconj))) +
  scale_x_continuous(limits = limitseigvalRe) +
  scale_y_continuous(limits = limitseigvalIm) +
  geom_point() +
  theme(legend.position = "bottom") +
  labs(x = "Real part dominant eigenvalue (with conjugation)",
       y = "Imaginary part dominant eigenvalue (with conjugation)",
       caption = paste(niter, "iterations")) +
  facet_grid(rows = vars(nspecies, conjratecode),
             cols = vars(taxmatcode, abunmodelcode, cost),
             labeller = mylabeller)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "eigenvaluesdistrconj.png"))
}

if(simulateinvasion == TRUE) {
  CreatePlot(fillvar = "infgrowthfrac", filltitle = "Fraction infinite\ngrowth",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedfrac", filltitle = "Fraction equilibrium\nreached",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "tmaxshortfrac", filltitle = "Fraction tmax too short",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "infgrowthconjfrac",
             filltitle = "Fraction infinite growth\nwith conjugation",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedconjfrac",
             filltitle = "Fraction equilibrium\nreached with\nconjugation",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "tmaxshortconjfrac",
             filltitle = "Fraction tmax too short\nwith conjugation",
             filltype = "continuous", limits = limitsfraction)
}
CreatePlot(fillvar = "fracreal", filltitle = "Fraction real",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracrealconj",
           filltitle = "Fraction real\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fraceigvalRep", filltitle = "Fraction repeated eigenvalues",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fraceigvalRepconj",
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
                                  near(conjratecode, 1), near(taxmatcode, 1))

ggplot(data = datatotalfiltercostconj, aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = selfintmean), size = 1) +
  facet_grid(species + nspecies ~ abunmodelcode, labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations"))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsintmean.png"))
}

ggplot(data = datatotalfiltercostconj, aes(x = selfintmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = intmean), size = 1) +
  facet_grid(species + nspecies ~ abunmodelcode, labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations"))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsselfintmean.png"))
}

# Calculate density if only intraspecies interactions are present
CreatePlot(dataplot = datatotalfiltercostconj, fillvar = "growthrate/selfintmean",
           filltitle = "Mean carrying capacity", filltype = "continuous",
           facetx = "abunmodelcode", facety = "nspecies",
           filename = "carrycap")
CreatePlot(dataplot = datatotalfiltercostconj, fillvar = "growthrate/selfintmean",
           filltitle = "Mean carrying capacity", filltype = "continuous",
           facetx = "abunmodelcode", facety = "nspecies",
           rotate_legend = TRUE, filename = "carrycaprotated")

## Plot summary data on the number of iterations in creating intmat needed to
# find a stable equilibrium with the model without plasmids
CreatePlot(fillvar = "iterintmatmin", filltitle =
             paste("Minimum number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(1, niterintmat))
CreatePlot(fillvar = "iterintmatmean", filltitle = 
             paste("Mean number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(1, niterintmat))
CreatePlot(fillvar = "iterintmatmedian", filltitle = 
             paste("Median number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(1, niterintmat))
CreatePlot(fillvar = "iterintmatmax", filltitle = 
             paste("Maximum number of\niterations to reach\nstable equilibrium"),
           filltype = "continuous", limits = c(1, niterintmat))

if(simulateinvasion == TRUE) {
  subplasmidfree <- "Perturbation with plasmid-free bacteria"
  subplasmidbearing <- "Perturbation with plasmid-bearing bacteria"
  
  limitstime <- range(plotdata[, "timefinalmedian"],
                      plotdata[, "timefinalconjmedian"])
  title <- "Time to reach equilibrium after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  CreatePlot(fillvar = "timefinalmedian", filltitle = "Median time",
             filltype = "continuous", limits = limitstime,
             title = title, subtitle = subplasmidfree, rotate_legend = TRUE)
  CreatePlot(fillvar = "timefinalconjmedian", filltitle = "Median time",
             filltype = "continuous", limits = limitstime,
             title = title, subtitle = subplasmidbearing, rotate_legend = TRUE)
  
  ## Plots on survival and extinction after perturbation
  limitsmeannspecies <- range(plotdata[, "npopRmean"],
                          plotdata[, "nspeciesconjmean"], finite = TRUE)
  title <- "Number of species surviving after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  # Note: In the model without plasmids, the number of populations is equal to
  # the number of species.
  CreatePlot(fillvar = "npopRmean", filltitle = "Mean total number\nof species",
             filltype = "continuous", limits = limitsmeannspecies,
             title = title, subtitle = subplasmidfree,
             filename = "nspeciesRmean")
  CreatePlot(fillvar = "nspeciesconjmean", filltitle = "Mean total number\nof species",
             filltype = "continuous", limits = limitsmeannspecies,
             title = title, subtitle = subplasmidbearing)
  
  title <- "Number of species going extinct after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                               near(conjratecode, 1), near(taxmatcode, 1)),
             fillvar = "nspecies - npopRmean",
             filltitle = "Mean number of\nspecies going extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidfree,
             facetx = "abunmodelcode", facety = "nspecies",
             filename = "nspeciesRmeanextinctfewfacets")
  CreatePlot(fillvar = "nspecies - npopRmean",
             filltitle = "Mean number of\nspecies going extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidfree,
             filename = "nspeciesRmeanextinct")
  CreatePlot(fillvar = "nspecies - nspeciesconjmean",
             filltitle = "Mean number of\nspecies going extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidbearing,
             filename = "nspeciesconjmeanextinct")
  
  title <- "Fraction of species going extinct after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                               near(conjratecode, 1), near(taxmatcode, 1)),
             fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of\nspecies going extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree,
             facetx = "abunmodelcode", facety = "nspecies",
             filename = "fracspeciesextinctmeanfewfacets")
  CreatePlot(fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of\nspecies going extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree,
             filename = "fracspeciesextinctmean")
  CreatePlot(fillvar = "1 - (nspeciesconjmean / nspecies)",
             filltitle = "Mean fraction of\nspecies going extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing,
             filename = "fracspeciesextinctmeanconj")
  
  title <- "Fraction of species after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  CreatePlot(fillvar = "fracspeciesRconjmean",
             filltitle = "Mean fraction of surviving species\nthat have a plasmid-free population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "npopRconjmean / nspecies",
             filltitle = "Mean fraction of initial species\nthat have a plasmid-free population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "fracspeciesPconjmean",
             filltitle = "Mean fraction of surviving species\nthat have a plasmid-bearing population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "npopPconjmean / nspecies",
             filltitle = "Mean fraction of initial species\nthat have a plasmid-bearing population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  
  ## Plots of fractions of bacteria that are plasmid-free or plasmid-bearing
  # after perturbations. Only abundances where equilibrium was reached are
  # considered.
  title <- "Fraction bacteria after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  CreatePlot(fillvar = "relabunRconjmin",
             filltitle = "Minimum fraction of bacteria\nthat is plasmid-free",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "relabunRconjmean",
             filltitle = "Mean fraction of bacteria\nthat is plasmid-free",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "relabunRconjmedian",
             filltitle = "Median fraction of bacteria\nthat is plasmid-free",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "relabunRconjmax",
             filltitle = "Maximum fraction of bacteria\nthat is plasmid-free",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  
  CreatePlot(fillvar = "relabunPconjmin",
             filltitle = "Minimum fraction of bacteria\nthat is plasmid-bearing",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "relabunPconjmean",
             filltitle = "Mean fraction of bacteria\nthat is plasmid-bearing",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "relabunPconjmedian",
             filltitle = "Median fraction of bacteria\nthat is plasmid-bearing",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "relabunPconjmax",
             filltitle = "Maximum fraction of bacteria\nthat is plasmid-bearing",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)

  CreatePlot(fillvar = "fracPformedbypertpopmin",
             filltitle = paste("Minimum fraction of plasmid-bearing bacteria",
                               "\nbelonging to the initially plasmid-bearing strain"),
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "fracPformedbypertpopmean",
             filltitle = paste("Mean fraction of plasmid-bearing bacteria",
                               "\nbelonging to the initially plasmid-bearing strain"),
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "fracPformedbypertpopmedian",
             filltitle = paste("Median fraction of plasmid-bearing bacteria",
                               "\nbelonging to the initially plasmid-bearing strain"),
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "fracPformedbypertpopmax",
             filltitle = paste("Maximum fraction of plasmid-bearing bacteria",
                               "\nbelonging to the initially plasmid-bearing strain"),
             filltype = "continuous", limits = limitsfraction)

  ## Plot of total abundances after perturbation with plasmid-free bacteria in
  # models without plasmids. Only abundances where equilibrium was reached are
  # considered. Although costs and conjugation rates do not influence the
  # outcome, they are included as facets to facilitate comparison with plots of
  # abundances after perturbation with plasmid-bearing bacteria.
  
  # log10(0) leads to -Inf and is displayed as a gray square in the plots,
  # leading to different numbers of gray squares between plots. It is also
  # confusing because gray squares can also indicate no data is available
  # because all simulations were discarded because equilibrium was not reached.
  # Therefor I plot log10(1 + x) instead of log10(x). I do not use the built-in
  # function log1p(x) because that uses natural logarithms.
  title <- "Total abundance after perturbation"
  if(PInMostAbun == FALSE) {
    title <- paste(title, "(PinLeastAbun)")
  }
  CreatePlot(fillvar = "log10(1 + abuntotalmin)",
             filltitle = "Log10(1 + Minimum of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidfree)
  CreatePlot(fillvar = "log10(1 + abuntotalmean)",
             filltitle = "Log10(1 + Mean of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidfree)
  CreatePlot(fillvar = "log10(1 + abuntotalmedian)",
             filltitle = "Log10(1 + Median of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidfree)
  CreatePlot(fillvar = "log10(1 + abuntotalmax)",
             filltitle = "Log10(1 + Maximum of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidfree)

  ## Plot of total abundances after perturbation with plasmid-bearing bacteria in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  CreatePlot(fillvar = "log10(1 + abuntotalconjmin)",
             filltitle = "Log10(1 + Minimum total\nabundance)", filltype = "continuous",
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abuntotalconjmean)",
             filltitle = "Log10(1 + Mean total\nabundance)", filltype = "continuous",
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abuntotalconjmedian)",
             filltitle = "Log10(1 + Median total\nabundance)", filltype = "continuous",
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abuntotalconjmax)",
             filltitle = "Log10(1 + Maximum total\nabundance)", filltype = "continuous",
             title = title, subtitle = subplasmidbearing)
  
  ## Plot total abundances of plasmid-free populations after perturbations in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  CreatePlot(fillvar = "log10(1 + abunRtotalconjmin)",
             filltitle = "Log10(1 + Minimum of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abunRtotalconjmean)",
             filltitle = "Log10(1 + Mean of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abunRtotalconjmedian)",
             filltitle = "Log10(1 + Median of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abunRtotalconjmax)",
             filltitle = "Log10(1 + Maximum of plasmid-\nfree bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)  
  
  ## Plot total abundances of plasmid-bearing populations after perturbations for
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  CreatePlot(fillvar = "log10(1 + abunPtotalconjmin)",
             filltitle = "Log10(1 + Minimum of plasmid-\nbearing bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abunPtotalconjmean)",
             filltitle = "Log10(1 + Mean of plasmid-\nbearing bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abunPtotalconjmedian)",
             filltitle = "Log10(1 + Median of plasmid-\nbearing bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "log10(1 + abunPtotalconjmax)",
             filltitle = "Log10(1 + Maximum of plasmid-\nbearing bacteria)",
             filltype = "continuous", title = title, subtitle = subplasmidbearing)
  
  ## Plots comparing species abundances after perturbation with plasmid-bearing
  # bacteria
  limits <- range(c(plotdata[, "relabunRsp1median"],
                    plotdata[, "relabunconjsp1median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp1median",
             filltitle = "Median rel. abundance sp1 after\nperturbation with R1",
             filltype = "continuous", limits = limits)
  CreatePlot(fillvar = "relabunconjsp1median",
             filltitle = "Median rel. abundance sp1 after\nperturbation with P1",
             filltype = "continuous", limits = limits)
  
  CreatePlot(fillvar = "relabunRsp1median",
             filltitle = "Median rel. abundance sp1 after\nperturbation with R1",
             filltype = "continuous", limits = limitsfraction,
             filename = "relabunRsp1mediancontinuouschangedlim")
  CreatePlot(fillvar = "relabunconjsp1median",
             filltitle = "Median rel. abundance sp1 after\nperturbation with P1",
             filltype = "continuous", limits = limitsfraction,
             filename = "relabunconjsp1mediancontinuouschangedlim")
  
  limits <- range(c(plotdata[, "relabunRsp2median"],
                    plotdata[, "relabunconjsp2median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp2median",
             filltitle = "Median rel. abundance sp2 after\nperturbation with R1",
             filltype = "continuous", limits = limits)
  CreatePlot(fillvar = "relabunconjsp2median",
             filltitle = "Median rel. abundance sp2 after\nperturbation with P1",
             filltype = "continuous", limits = limits)
  
  limits <- range(c(plotdata[, "relabunRsp3median"],
                    plotdata[, "relabunconjsp3median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp3median",
             filltitle = "Median rel. abundance sp3 after\nperturbation with R1",
             filltype = "continuous", limits = limits)
  CreatePlot(fillvar = "relabunconjsp3median",
             filltitle = "Median rel. abundance sp3 after\nperturbation with P1",
             filltype = "continuous", limits = limits)
  
  limits <- range(c(plotdata[, "relabunRsp4median"],
                    plotdata[, "relabunconjsp4median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp4median",
             filltitle = "Median rel. abundance sp4 after\nperturbation with R1",
             filltype = "continuous", limits = limits)
  CreatePlot(fillvar = "relabunconjsp4median",
             filltitle = "Median rel. abundance sp4 after\nperturbation with P1",
             filltype = "continuous", limits = limits)
  
  CreatePlot(fillvar = "log10(1 + relabunRsp4median)",
             filltitle = "Log10(1 + Median rel. abundance sp4 after\nperturbation with R1)",
             filltype = "continuous")
  CreatePlot(fillvar = "log10(1 + relabunconjsp4median)",
             filltitle = "Log10(1 + Median rel. abundance sp4 after\nperturbation with P1)",
             filltype = "continuous")
  
  limits <- range(c(plotdata[, "relabunRsp5median"],
                    plotdata[, "relabunconjsp5median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp5median",
             filltitle = "Median rel. abundance sp5 after\nperturbation with R1",
             filltype = "continuous", limits = limits)
  CreatePlot(fillvar = "relabunconjsp5median",
             filltitle = "Median rel. abundance sp5 after\nperturbation with P1",
             filltype = "continuous", limits = limits)
  
  limits <- range(c(plotdata[, "relabunRsp6median"],
                    plotdata[, "relabunconjsp6median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp6median",
             filltitle = "Median rel. abundance sp6 after\nperturbation with R1",
             filltype = "continuous", limits = limits)
  CreatePlot(fillvar = "relabunconjsp6median",
             filltitle = "Median rel. abundance sp6 after\nperturbation with P1",
             filltype = "continuous", limits = limits)
  
  CreatePlot(fillvar = "log10(1 + relabunRsp6median)",
             filltitle = "Log10(1 + Median rel. abundance sp6 after\nperturbation with R1)",
             filltype = "continuous")
  CreatePlot(fillvar = "log10(1 + relabunconjsp6median)",
             filltitle = "Log10(1 + Median rel. abundance sp6 after\nperturbation with P1)",
             filltype = "continuous")
}

## Compare abundance models ##
compareabun <- NULL
for(abunmodel in abunmodelset) {
  for(nspecies in nspeciesset) {
    comparetemp <- data.frame(
      nspecies = as.factor(nspecies),
      species = as.factor(1:nspecies),
      abun = switch(abunmodel,
                    "brokenstick" = brokenstick(nspecies = nspecies,
                                                totalabun = totalabun,
                                                takelimit = TRUE),
                    "dompreempt" = dompreempt(nspecies = nspecies,
                                              totalabun = totalabun,
                                              takelimit = TRUE)),
      model = abunmodel)
    compareabun <- rbind(compareabun, comparetemp)
  }
}
compareabun[, "group"] <- paste0(compareabun[, "nspecies"],
                                 " species, ", compareabun[, "model"])
if(!exists("DateTimeStamp")) {
  DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
}

if(length(abunmodelset) > 1) {
  facet_rows <- vars(model)
} else {
  facet_rows <- NULL
}
plotcompareabun <- ggplot(data = compareabun,
                          aes(x = species, y = abun,
                              color = nspecies, lty = model)) +
  theme_bw(base_size = 14) +
  scale_x_discrete(limits = factor(1:maxnspecies)) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_grid(rows = facet_rows, cols = vars(nspecies), labeller = "label_both") +
  theme(panel.grid = element_line(size = 1),
        legend.position = "bottom", legend.just = "left",
        legend.margin = margin(-5, 0, -5, 0),
        legend.box = "vertical", legend.box.just = "left",
        legend.box.margin = margin(-5, 0, 0, 0)) +
  labs(title = "Comparing abundance models", subtitle = "Linear scale",
       x = "Species rank", y = "Species abundance") +
  geom_line(aes(group = group), size = 1.25) +
  geom_point(size = 3)
print(plotcompareabun)
if(saveplots == TRUE) {
  filename <- paste0(DateTimeStamp, "compareabunmodel.png")
  ggsave(filename, width = 16, height = 16, units = "cm")
}

plotcompareabunlog <- ggplot(data = compareabun,
                             aes(x = species, y = abun,
                                 color = nspecies, lty = model)) +
  theme_bw(base_size = 14) +
  scale_x_discrete(limits = factor(1:maxnspecies)) +
  scale_y_continuous(trans = "log10") +
  facet_grid(rows = facet_rows, cols = vars(nspecies), labeller = "label_both") +
  theme(panel.grid = element_line(size = 1),
        legend.position = "bottom", legend.just = "left",
        legend.margin = margin(-5, 0, -5, 0),
        legend.box = "vertical", legend.box.just = "left",
        legend.box.margin = margin(-5, 0, 0, 0)) +
  labs(title = "Comparing abundance models", subtitle = "Logarithmic scale",
       x = "Species rank", y = "Species abundance") +
  geom_line(aes(group = group), size = 1.25) +
  geom_point(size = 3)
print(plotcompareabunlog)
if(saveplots == TRUE) {
  filename <- paste0(DateTimeStamp, "compareabunmodellog.png")
  ggsave(filename, width = 16, height = 16, units = "cm")
}

write.csv(compareabun,
          file = paste0(DateTimeStamp, "compareabunmodel.csv"),
          quote = FALSE, row.names = FALSE)
