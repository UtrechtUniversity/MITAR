################################################################################
## Factors influencing the invasion of plasmid-free and plasmid-bearing       ##
## bacteria in a microbiome                                                   ##
################################################################################

# Test Git push

#### Introduction ####
# Use simulations with a generalised Lotka-Volterra model elaborated with
# plasmid-bearing populations and conjugation to investigate how the
# number of species and the evenness of their abundances in a microbiome
# influence the ability of bacteria to invade that microbiome, and to
# investigate how is this changed if the invading bacterium carries a plasmid,
# a (plasmid-free population) of the invading species is already present in the
# microbiome, and conjugation is, or is not, relatedness-dependent.

# NOTE:
# - Some functions differ from the functions for the case where the plasmid is
#   introduced through an already-present species.

#### References ####
# Alderliesten JB, Duxbury SJN, Zwart MP, de Visser JAGM, Stegeman S, Fische
# EAJ. 2020. Effect of donor-recipient relatedness on the plasmid conjugation
# frequency: a meta-analysis. BMC Microbiology 20:135.

# Edelstein-Keshet L. 2005. Mathematical models in biology. Society for
# industrial and applied mathematics.

# Lischke H, Löffler TJ. 2017. Finding all multiple stable fixpoints of n-species
# Lotka-Volterra competition models. Theoretical Population Biology 115:24-34.

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
# into account. A check for repeated eigenvalues is included, but note that the
# pair a +/- bi are not repeated eigenvalues. See p. 133 of Edelstein-Keshet
# 2005 and section 10.4.3.3 from paste0("https://eng.libretexts.org/",
# "Bookshelves/Industrial_and_Systems_Engineering/Book%3A_Chemical_Process",
# "_Dynamics_and_Controls_(Woolf)")


#### Optionally to do ####
# Supplying the analytic Jacobian to the solver in perturbequilibrium() might
# speed up the integration. See Box 1 in Lischke 2017 for the structure of the
# Jacobian.

# The .groups argument in dplyr::summarise() is experimental, so maybe should be
# dropped. However, excluding it leads to messages every time it is called.

# To check if enough iterations are used: run several times for niter iterations
# and use more iterations if the variation in outcomes is too large.

# The matrix 'data' is converted to a tibble to efficiently get summary
# statistics. The tibble with summary statistics is then converted to a matrix
# to fill in part of the matrix plotdata. Using a tibble for 'data' from the
# start prevents some of this type-conversion, might make naming columns on
# assignment clearer, and enable use of the more efficient dplyr::bind_cols()
# instead of base::cbind(). However, originally using a data.frame instead of
# matrix to store data made progress really slow, so I should check how timing
# is affected when 'data', 'plotdata' or both are tibbles from the start.


#### Loading required packages ####
library(deSolve)   # checkequilibrium() and perturbequilibrium() call ode() if
                   # showplot and simulateinvasion are 'TRUE', respectively
library(dplyr)     # across(), full_join(), group_by(), near(), summarise()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()
# On the pipe operator (%>%), see ?'%>%' and Ch. 18 'Pipes' in Wickham 2017.


#### Settings and defining parameter space ####
# If simulateinvasion == TRUE, simulations over time are performed.
# If states become smaller than smallstate during the integration, they are set
# to 0.
# If the sum of absolute rates of change is equal to smallchange, equilibrium is
#   assumed to be reached and the integration is terminated.
# See the functions that use the arguments for more detailed info.

# Note: to simulate that each species belongs to a different class, an
# additional taxmattype has to be added. Then only the interspecies conjugation
# rates should be reduced, as the intraspecies conjugation rates of all species
# are still those given in conjrateset. This unchanged intraspecies conjugation
# rate makes it different from reducing conjrateset 1000-fold and using
# 'taxmatsame' as taxmattype.

## Basis parameter set
bifurparms <- FALSE
saveplots <- TRUE
niterintmat <- 1
smallstate <- 1e-3
finalsmallstate <- 1
smallchange <- 1e-2
totalabun <- 1e11
nspeciesset <- c(2, 4, 8, 16)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
costset <- c(0.05, 0.09)
costtype <- "absolute"
conjrateset <- list(rep(1e-12, maxnspecies))
# If taxmattype is "SameSpecies", the conjugation rate is the same for all
# populations. If taxmattype is "OtherClass", the interspecies conjugation rate
# to and from the initially plasmid-bearing population is reduced a 1000-fold.
taxmattypeset <- c("SameSpecies", "OtherClass")
# If PReplMostAbun is TRUE, some plasmid-free bacteria of the most-abundant
# species (i.e., species 1) are added to that species to simulate perturbation
# by plasmid-free bacteria, and some plasmid-free bacteria of the most-abundant
# species are replaced with plasmid-bearing bacteria of that species to simulate
# perturbation by plasmid-bearing bacteria.
# If PReplMostAbun is FALSE, some plasmid-free bacteria of the least-abundant
# species (i.e., species nspecies) are added to that species to simulate
# perturbation by plasmid-free bacteria, and some plasmid-free bacteria of the
# least-abundant species are replaced with plasmid-bearing bacteria of that
# species to simulate perturbation by plasmid-bearing bacteria.
PReplMostAbun <- TRUE
# To plot 16 species need 16 colours, currently only 11 so repeat them.
mycol <- rep(c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
               "darkorange", "green1", "yellow", "hotpink"), 2)

## Parameters for detailed local stability analysis, not simulating invasion
niter <- 100
simulateinvasion <- FALSE
intmeanset <- seq(from = -1e-11, to = 5e-12, by = 1e-12)
selfintmeanset <- seq(from = -1e-11, to = 0, by = 1e-12)

## Smaller parameter set to simulate invasion
niter <- 50
simulateinvasion <- TRUE
intmeanset <- seq(from = -1e-11, to = 5e-12, by = 1e-12)
selfintmeanset <- seq(from = -1e-11, to = 0, by = 1e-12)

## Small parameter set to test code
nspeciesset <- c(2, 8)
maxnspecies <- max(nspeciesset)
conjrateset <- list(rep(1e-13, maxnspecies), rep(1e-12, maxnspecies))
niter <- 2
simulateinvasion <- TRUE
intmeanset <- c(-6e-12, 6e-12)
selfintmeanset <- c(-1.2e-11, -6e-12)

## Larger test set
nspeciesset <- c(2, 8)
maxnspecies <- max(nspeciesset)
conjrateset <- list(rep(1e-13, maxnspecies), rep(1e-12, maxnspecies))
niter <- 25
simulateinvasion <- TRUE
intmeanset <- c(1e-11, -6e-12, 0, 6e-12, 1e11)
selfintmeanset <- c(-1e-11, -6e-12, -1e-12)

## Parameter set to create bifurcation-like plots showing the border of
# epidemiological stability in the conjugation rate/cost space
bifurparms <- TRUE
nspeciesset <- c(2, 16) # Check if 4 or 8 species is the same as 16 species.
maxnspecies <- max(nspeciesset)
costset <- seq(from = 0, to = 0.2, by = 0.001)
seqconjrate <- 10^seq(from = -13.5, to = -11.5, by = 0.01)
conjrateset <- NULL
# The conjugation rate given here is the 'overall' conjugation rate. For all
# species, the intraspecies conjugation rate will be equal to this overall
# conjugation rate. The interspecies conjugation rate will be equal to this
# overall conjugation rate if all populations belong to the same species.
# If the initially plasmid-bearing species is considered to be from another
# taxonomic class than the initially plasmid-free species, these overall
# conjugation rates are reduced a thousand-fold to obtain the interspecies
# conjugation rate between initially plasmid-bearing and initially plasmid-free
# species.
for(conjrate in seqconjrate) {
  conjrateset <- c(conjrateset, list(rep(conjrate, maxnspecies)))
}
niter <- 1
simulateinvasion <- FALSE
intmeanset <- 0 # c(-1e-11, 0, 1e-11)
selfintmeanset <- -0.5e-11 # c(-1e-11, -0.5e-11, 0)


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
    S0 <- n[seq_len(nspecies)]
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
    for(iterindex in seq_len(niterabun)) {
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

brokenstick_fast <- function(nspecies, totalabun, takelimit = NULL,
                             niterabun = 1000) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
            length(totalabun) == 1, totalabun > 0)
  abunmat <- matrix(data = NA, nrow = niterabun, ncol = nspecies)
  for(iterindex in seq_len(niterabun)) {
    # sorting to get positive differences in the next step
    breakpoints <- sort(runif(nspecies - 1, min = 0, max = totalabun))
    # sorting because otherwise taking the mean does not make sense
    abunmat[iterindex, ] <- sort(diff(c(0, breakpoints, totalabun)))
  }
  sort(colMeans(abunmat), decreasing = TRUE)
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
    abun <- totalabun*0.75*0.25^(seq_len(nspecies) - 1)
  } else {
    remainingabun <- totalabun
    for(speciesindex in seq_len(nspecies)) {
      abun.temp <- runif(1, min = 0.5, max = 1) * remainingabun
      abun[speciesindex] <- abun.temp
      remainingabun <- remainingabun - abun.temp
    }
  }
  # Scaling abundances to get specified total abundance
  abun <- abun / (1 - (0.25^nspecies))
  return(abun)
}

dompreempt_fast <- function(nspecies, totalabun, takelimit = NULL) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
            length(totalabun) == 1, totalabun > 0)
  abun <- rep(NA, nspecies)
  abun <- totalabun*0.75*0.25^(seq_len(nspecies) - 1)
  
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
# drawn from the distribution 'intdistr'. Diagonal entries are intraspecies
# interactions drawn from the distribution 'selfintdistr', which is truncated to
# obtain negative intraspecies interactions. To get fixed values for the
# interactions, choose the uniform distribution and provide the desired value
# both as the minimum and maximum of the range, or provide the standard
# deviation of the normal distribution as zero. The other arguments specify the
# distributions from which interaction coefficients are drawn. See getintmat()
# for a faster, more succinct version of this function.
getintmat_elaborate <- function(nspecies, sparsity = 0,
                                intdistr = "normal", intmean = 0, intsd = 5e-12,
                                intrange = c(-1e-11, 1e-11),
                                selfintdistr = "normal",
                                selfintmean = -5e-12, selfintsd = 9e-12,
                                selfintrange = c(-1e-11, 0)) {
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
           intcoefs <- NULL
           stop("'intdistr' should be 'normal' or 'uniform'.")
         }
  )
  intmat <- matrix(intcoefs, nrow = nspecies, ncol = nspecies)
  
  switch(selfintdistr,
         normal = {
           stopifnot(length(selfintmean) == 1, length(selfintsd) == 1,
                     selfintsd >= 0)
           # Draw variates from a truncated normal distribution to ensure
           # negative intraspecies interactions without having to redraw for
           # positive deviates, using rtnorm() from the package TruncatedNormal.
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
           intmat <- NULL
           stop("'selfintdistr' should be 'normal' or 'uniform'.")
         }
  )
  
  if(sparsity > 0) {
    nsparseint <- round(sparsity*(nspecies^2 - nspecies))
    # Create indexmatrix with column- and row indices
    indexmat <- matrix(c(rep(seq_len(nspecies), nspecies),
                         rep(seq_len(nspecies), each = nspecies)),
                       ncol = 2, dimnames = list(NULL, c("row", "column")))
    # Only keep rows of indexmat specifying off-diagonal entries of intmat, to
    # ensure that intraspecies interaction coefficients never become sparse.
    indexmat <- indexmat[indexmat[, "row"] != indexmat[, "column"], , drop = FALSE]
    # Sample rows of indexmat to get index of matrix entries that become sparse
    sparse_index <- sample(seq_len(dim(indexmat)[1]), nsparseint)
    # Only keep rows of indexmat that were drawn from the sample. Use drop =
    # FALSE to prevent sparse_index becoming a vector of length two such that
    # two elements are set to zero if nsparseint == 1. 
    indexmat <- indexmat[sparse_index, , drop = FALSE]
    # Set entries of the interaction matrix indicated by indexmat to zero.
    intmat[indexmat] <- 0
  }
  return(intmat)
}

# Create a matrix of scaled interaction coefficients for nspecies species with
# units mL cell^-1 h^-1, where element aij represents cij * ri / Ki in the
# textbook-notation (see comments on the generalised Lotka-Volterra model given
# above), such that elements aii on the diagonal equal ri / Ki, with Ki being
# the carrying capacity of species i. Off-diagonal entries are interspecies
# interaction coefficients drawn from a normal distribution with mean and
# standard deviation given by intmean and intsd, respectively. Diagonal entries
# are intraspecies interactions drawn from a truncated normal distribution with
# mean and standard deviation given by selfintmean and selfintsd, respectively,
# and is truncated to obtain negative intraspecies interactions. See function
# 'getintmat_elaborate()' for a more elaborate but slower version of this
# function.
getintmat <- function(nspecies, intmean = 0, intsd = 5e-12,
                      selfintmean = -5e-12, selfintsd = 9e-12) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
            length(intmean) == 1, length(intsd) == 1, intsd >= 0,
            length(selfintmean) == 1, length(selfintsd) == 1, selfintsd >= 0)
  
  intmat <- matrix(rnorm(n = nspecies^2, mean = intmean, sd = intsd),
                   nrow = nspecies, ncol = nspecies)
  
  # Draw variates from a truncated normal distribution to ensure
  # negative intraspecies interactions without having to redraw for
  # positive deviates, using rtnorm() from the package TruncatedNormal.
  diag(intmat) <- rtnorm(n = nspecies, mu = selfintmean, sd = selfintsd,
                         lb = -Inf, ub = 0, method = "invtransfo")
  return(intmat)
}

# Calculate the required growth rates to obtain the specified species abundances
# at equilibrium given the interaction matrix. The growth rate can be negative
# if growth is assumed to be slower than washout. See getgrowthrate() for a
# faster, more succinct version of this function.
getgrowthrate_elaborate <- function(abundance, intmat) {
  stopifnot(length(abundance) == dim(intmat)[1])
  growthrate <- -intmat %*% abundance
  if(any(is.complex(growthrate))) {
    countcomplex <- length(which(Im(growthrate) != 0))
    warntext <- paste(countcomplex, "growth rates contain an imaginary part.")
    warning(warntext)
  }
  return(growthrate)
}

# Calculate the required growth rates to obtain the specified species abundances
# at equilibrium given the interaction matrix. The growth rate can be negative
# if growth is assumed to be slower than washout. See getgrowthrate_elaborate()
# for a more elaborate version of this function.
getgrowthrate <- function(abundance, intmat) {
  stopifnot(length(abundance) == dim(intmat)[1])
  -intmat %*% abundance
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
    out <- ode(y = abundance, times = times, func = gLV,
               parms = list(growthrate = growthrate, intmat = intmat),
               rootfunc = rootfun,
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
#     and, if taxmattypeset contains "OtherClass", "OtherClass" are used.
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
#     Only the submatrix taxmat[seq_len(nspecies), seq_len(nspecies)] is used if
#     nspecies < max(nspeciesset). The use of multiple values for taxmat is
#     NOT supported.
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
  conjratensp <- conjrate[seq_len(nspecies)]
  taxmatnsp <- taxmat[seq_len(nspecies), seq_len(nspecies)]
  
  # To obtain interspecies conjugation rates for the different levels of
  # taxonomic relatedness between donor and recipients, the intraspecies
  # conjugation rates are multiplied with the following conversion factors
  # (based on the decrease in conjugation frequencies described in Alderliesten
  # 2020).
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

# Can be used if taxmat contains nothing else than 'SameSpecies' or 'OtherClass'
getconjmat_fast <- function(nspecies, conjrate, taxmat) {
  stopifnot(all(diag(taxmat) == "SameSpecies"),
            isSymmetric.matrix(unname(taxmat)))
  sel_sp <- seq_len(nspecies)

  # To obtain interspecies conjugation rates for the species belonging to
  # another taxonomic class than the donor, the intraspecies conjugation rates
  # are multiplied with 0.001 (based on the decrease in conjugation frequencies
  # described in Alderliesten 2020).
  convmat <- matrix(1, nrow = nspecies, ncol = nspecies)
  convmat[which(taxmat[sel_sp, sel_sp] == "OtherClass")] <- 0.001
  
  # Multiply column n of convmat giving the conversion factors for conjugation
  # from donor species n to the different recipient species, with element n of
  # conjrate.
  convmat %*% diag(conjrate[sel_sp])
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
           jaclow <- jac[(nspecies + 1):(2*nspecies), seq_len(nspecies)]
           if(!isTRUE(all.equal(range(jaclow), c(0, 0), check.attributes = FALSE))) {
             warning("Jacobian matrix does not contain a block of zeros in the",
                     " lower-left corner,\nso currently used determination of",
                     " ecological stability is invalid")
             print(jaclow)
           }
           eigval <- eigen(x = jac[seq_len(nspecies), seq_len(nspecies)],
                           symmetric = FALSE, only.values = TRUE)$values
         },
         "epi" = {
           jac <- jacobian.full(y = abundance, func = gLVConj,
                                parms = list(growthrate = growthrate, intmat = intmat,
                                             cost = cost, conjmat = conjmat))
           indexP <- (nspecies + 1):(2*nspecies)
           jaclow <- jac[indexP, seq_len(nspecies)]
           if(!isTRUE(all.equal(range(jaclow), c(0, 0), check.attributes = FALSE))) {
             warning("Jacobian matrix does not contain a block of zeros in the",
                     " lower-left corner,\nso currently used determination of",
                     " epidemiological stability is invalid")
             print(jaclow)
           }
           eigval <- eigen(x = jac[indexP, indexP],
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
#   vector. If the latter are not 0, a warning is issued.
# intmat, growthrate, and conjmat are supposed to be the output of getintmat(),
#   getgrowthrate(), and getconjmat()
# cost is a vector of plasmid costs in growth rate. 
# model should be 'gLV' (no plasmids modelled) or 'gLVConj' (to include plasmid-
#   bearing populations in the model).
# pertpop is a character vector with the name(s) of the population(s) to be
#   perturbed, e.g., pertpop = "R1" or pertpop = c("R1", "P1").
# pertmagn gives the absolute increase in populations for the perturbation
# tmax and tstep determine the maximum timestep and the intervals at which
#   abundances should be calculated. Timesteps of 1 will be used until the first
#   timestep, to reduce the chance using too large timesteps at the start of the
#   integration. Since ode() uses variable step-size methods, those are not the
#   only times that integration occurs.
# If showplot == TRUE, the result is plotted (which is slow) or (if saveplots
#   == TRUE) saved to a .png file.
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
#   very large. A warning is issued (unless silentinfgrowth is FALSE) and the
#   variable 'infgrowth' is set to 1 if this occurs.
# When the simulation is finished, abundances smaller than finalsmallstate with
#   a negative derivative (i.e., small abundances that are declining) are set to
#   zero
# Output:
# Returns abunfinal, a list containing:
#   - R with absolute abundances of plasmid-free populations for each species,
#   - Rtotal with the sum of abundances of plasmid-free populations
#   - npopR with the number of plasmid-free populations
#   - P with absolute abundances of plasmid-bearing population for each species
#     (or NULL if model == "gLV")
#   - Ptotal with the sum of abundances of plasmid-bearing populations (or NULL
#     if model == "gLV")
#   - npopP with the number of plasmid-bearing populations
#   - pertpopconjsurvived indicating the initially plasmid-bearing population
#     survived.
#   - timepertpopconjextinct indicating the last recorded time the initially
#     plasmid-bearing population was extant, i.e., abundance > finalsmallstate
#   - timefinal indicating the last recorded time
#   - tmaxshort indicating tmax was not reached but no infinite growth occured
#   - eqreached indicating if equilibrium has been reached
#   - infgrowth indicating if infinite growth was detected
# Notes
# - Although this function returns absolute abundances, various of these are
#   converted to relative abundances later on in the script.
# - The current implementation of pertpopconjsurvived does NOT handle
#   perturbation by multiple populations.
# - See perturbequilibrium_gLV_fast() and perturbequilibrium_gLVConj_fast for
#   faster implementations for those models.
perturbequilibrium <- function(abundance, intmat, growthrate, cost, conjmat,
                               model, pertpop, pertmagn = 1000,
                               tmax = 5e3, tstep = 10, showplot = TRUE,
                               plotepistabwithP = FALSE, verbose = FALSE,
                               silentinfgrowth = FALSE,
                               silenteqnotreached = FALSE) {
  
  # Name abundances, set line type and colors, get derivatives of initial state
  # in the plasmid-free model.
  if(model == "gLV") {
    nspecies <- length(abundance)
    index_sp <- seq_len(nspecies)
    indexR <- seq_len(nspecies)
    names(abundance) <- paste0(rep("R", nspecies), index_sp)
    lty <- 1
    col <- mycol[index_sp]
    derivatives <- unlist(
      gLV(t = 0, n = abundance,
          parms = list(growthrate = growthrate, intmat = intmat)
      )
    )
  }
  
  if(model == "gLVConj") {
    nspecies <- length(abundance)/2
    index_sp <- seq_len(nspecies)
    indexR <- seq_len(nspecies)
    indexP <- (nspecies + 1):(2*nspecies)
    names(abundance) <- c(paste0(rep("R", nspecies), index_sp),
                          paste0(rep("P", nspecies), index_sp))
    lty <- rep(c(1, 2), each = nspecies)
    col <- rep(mycol[index_sp], 2)
    if(!all(near(abundance[indexP], 0))) {
      warning("Initial state is NOT plasmid-free.")
    }
    # Derivatives of model without plasmids, to check if initial state is an
    # equilibrium in that model
    derivatives <- unlist(
      gLV(t = 0, n = abundance[indexR],
          parms = list(growthrate = growthrate, intmat = intmat)
      )
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
  
  # Create perturbed abundances. In the model with plasmid-bearing bacteria, the
  # equilibrium is perturbed by replacing plasmid-free bacteria with
  # plasmid-bearing bacteria, instead of adding plasmid-bearing bacteria to the
  # plasmid-free equilibrium. In that way the number of bacteria for each
  # species remains the same, such that the ecological equilibrium is only
  # affected by plasmid costs altering growth rates, not by changes in the
  # abundance of species. In the model with only plasmid-free bacteria (i.e.,
  # 'gLV'), this does not work because there are no plasmid-bearing bacteria to
  # replace plasmid-free bacteria, so with that model (additional) plasmid-free
  # bacteria are added to the plasmid-free equilibrium.
  abuninit <- abundance
  abunpert <- abundance
  abunpert[pertpop] <- abunpert[pertpop] + pertmagn
  
  if(model == "gLVConj") {
    # Some plasmid-free bacteria of the most-abundant (least-abundant) species
    # are replaced with plasmid-bearing bacteria of that species to simulate
    # perturbation by plasmid-bearing bacteria if PReplMostAbun is TRUE (FALSE).
    if(PReplMostAbun == TRUE) {
      pertpopminus <- "R1"
    } else {
      pertpopminus <- paste0("R", nspecies)
    }
    abunpert[pertpopminus] <-  abunpert[pertpopminus] - pertmagn
    
    if(any(abunpert < 0)) {
      warning("Initial abundances of some populations would become negative",
              " because\nthe pertubation is larger than the number of bacteria",
              " initially present.\nThese abundances have been set to zero",
              " instead.")
      abunpert[which(abunpert < 0)] <- 0
    }
  }
  
  if(verbose == TRUE) {
    print("abunpert =")
    print(abunpert)
  }
  
  # Perturb equilibrium
  times <- c(0:(tstep - 1), seq(from = tstep, to = tmax, by = tstep))
  if(model == "gLV") {
    out <- ode(y = abunpert, times = times, func = gLV,
               parms = list(growthrate = growthrate, intmat = intmat),
               rootfunc = rootfun,
               events = list(func = eventfun, root = TRUE, terminalroot = c(1, 2)))
  }
  if(model == "gLVConj") {
    out <- ode(y = abunpert, times = times, func = gLVConj,
               parms = list(growthrate = growthrate, intmat = intmat,
                            cost = cost, conjmat = conjmat),
               rootfunc = rootfunconj,
               events = list(func = eventfun, root = TRUE, terminalroot = c(1, 2)))
  }
  final <- tail(out, 1)
  timefinal <- final[, "time"]
  abunfinaltemp <- final[, -which(colnames(final) == "time")]
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
                      pertpopconjsurvived = NULL,
                      timepertpopconjextinct = NULL,
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
      pertpopconjsurvived <- NA
      timepertpopconjextinct <- NA
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
      
      if(abunfinaltemp[pertpop] > finalsmallstate) {
        pertpopconjsurvived <- 1
        timepertpopconjextinct <- NA
      } else {
        pertpopconjsurvived <- 0
        timepertpopconjextinct <- unname(out[max(
          which(out[, pertpopconj] >= finalsmallstate)), "time"])
      }
    }
    
    abunfinal <- list(R = abunfinaltemp[indexR],
                      Rtotal = sum(abunfinaltemp[indexR]),
                      npopR = npopR,
                      P = abunfinaltemp[indexP],
                      Ptotal = sum(abunfinaltemp[indexP]),
                      npopP = npopP, pertpopconjsurvived = pertpopconjsurvived,
                      timepertpopconjextinct = timepertpopconjextinct,
                      timefinal = timefinal, tmaxshort = tmaxshort,
                      eqreached = eqreached, infgrowth = infgrowth)
  }
  
  # eqinfoepi is created outside this function before this function is called.
  if(showplot == TRUE ||
     (plotepistabwithP == TRUE && model == "gLVConj" &&
      !is.na(abunfinal["Ptotal"]) && abunfinal["Ptotal"] > finalsmallstate &&
      eqinfoepi["eigvalRe"] < 0)) {
    subtitle <- paste0(abunmodel, ", intmean: ", intmean,
                       ", selfintmean: ", selfintmean, "\ncost: ", cost,
                       ", conjratecode: ", conjratecode,
                       ", taxmatcode: ", taxmatcode, ", iter: ", iter)
    if(saveplots == TRUE) {
      filename <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M_%OS3"), ".png")
      png(filename = filename, width = 480, height = 785)
    }
    matplot.deSolve(out, lty = lty, col = col, ylab = "Abundance",
                    log = "y", lwd = 2, legend = list(x = "bottomright"))
    grid()
    abline(h = abuninit)
    mtext(side = 1, line = -1, at = 0, adj = 0, cex = 0.9, subtitle)
    if(saveplots == TRUE) {
      dev.off()
      message("Saved plot ", filename)
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
getfracnotzero <- list(
  frac = ~if(any(!is.na(.x))) {mean(.x != 0, na.rm = TRUE)} else {NA}
)

# Function to create plots. The plotted object is returned, such that it can be
# further modified like any other ggplot object.
CreatePlot <- function(dataplot = plotdata, xvar = "intmean", yvar = "selfintmean",
                       contour_var = NULL, contour_breaks = 0.5, contour_col = NULL,
                       contour_lty = NULL,
                       limits = NULL, limx = NULL, limy = NULL, ratio = 1,
                       fillvar, filltitle, filltype = "discrete",
                       filllabels = NULL,
                       title = NULL, subtitle = NULL,
                       labx = "Mean interspecies interaction coefficient",
                       laby = "Mean intraspecies interaction coefficient",
                       tag = NULL, addstamp = FALSE, diagonal = "none",
                       linezero = TRUE,
                       facetx = "abunmodelcode + cost + taxmatcode",
                       facety = "conjratecode + nspecies",
                       dropfacets = TRUE,
                       as.table = TRUE,
                       marginx = NULL, marginy = NULL, base_size = 11,
                       rotate_x_labels = TRUE, rotate_legend = FALSE,
                       save = saveplots, width = 1650, height = 2675,
                       filename = NULL) {
  caption <- paste(unique(dataplot$niter), "iterations")
  if(exists("DateTimeStamp") == FALSE) {
    DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
    if(addstamp == TRUE) {
      warning("DateTimeStamp created to include in plot does not correspond to",
              " the filename of the dataset.")
    }
  }
  if(addstamp == TRUE) {
    caption <- paste0(caption, ", ", DateTimeStamp)
  }
  
  # Facets with only a single unique value are not shown if dropfacets = TRUE
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
    p <- p + geom_abline(intercept = 0, slope = -1, col = "grey", size = 1.1)
  }
  if(diagonal == "both" || diagonal == "minor") {
    p <- p + geom_abline(intercept = 0, slope = 1, col = "grey", size = 1.1)
  }
  if(linezero == TRUE) {
    p <- p + geom_vline(xintercept = 0, col = "grey", size = 1.1)
  }
  if(save == TRUE) {
    if(is.null(filename)) {
      filename <- paste0(fillvar, filltype)
    }
    
    # Add DateTimeStamp, remove spaces, replace non-alphanumeric characters with
    # underscores, and add the extension to create a valid file name.
    filename <- paste0(DateTimeStamp, filename)
    filename <- gsub(" ", "", filename)
    filename <- gsub("[^[:alnum:]_]", "_", filename)
    filename <- paste0(filename, ".png")
    
    if(file.exists(filename)) {
      warning("File '", filename, "' already exists, plot is not saved again!",
              call. = FALSE)
    } else {
      ggsave(filename, width = width, height = height, units = "px", dpi = 300)
    }
  }
  return(p)
}


#### Testing functions ####
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
# 
# # geteqinfo returns eigvalRe, eigvalIm, eigvalReSign, eigvalImSign, and eigvalRep
# geteqinfo(model = "gLV", abundance = abunbrokenstick, intmat = intmat1,
#           growthrate = growthratebrokenstick)
# geteqinfo(model = "gLV", abundance = abundompreempt, intmat = intmat1,
#           growthrate = growthratedompreempt)
# taxmat <- matrix(rep("SameSpecies", nspecies^2), nrow = nspecies,
#                  ncol = nspecies, byrow = TRUE)
# conjmat <- getconjmat(nspecies = nspecies, conjrate = rep(1e-12, nspecies),
#                       taxmat = taxmat)


#### Running the simulations ####
set.seed(seed = 314, kind = "default", normal.kind = "default", sample.kind = "default")
starttime <- Sys.time()

# Create matrix to store data
# To do:
# - Create variable 'newgrowthratecode' as NA and let it be written to the data
#   such that the same number of columns are present for this and for PinNew.
nrowplotdata <- prod(lengths(list(nspeciesset, abunmodelset, intmeanset,
                                  selfintmeanset, costset, conjrateset, taxmattypeset),
                             use.names = FALSE))
ncolplotdata <- if(simulateinvasion == TRUE) {
  18*4 + 4*4*maxnspecies + 38
} else {
  3*4 + 29
}
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = ncolplotdata)
nrowdatatotal <- prod(lengths(list(abunmodelset,intmeanset, selfintmeanset,
                                   costset, conjrateset, taxmattypeset),
                              use.names = FALSE))*niter*sum(nspeciesset)
datatotal <- matrix(data = NA, nrow = nrowdatatotal, ncol = 59 + 4*maxnspecies)
indexdatatotal <- 1

# Run simulations
rowindexplotdata <- 1
rowindexdata <- 1

for(nspecies in nspeciesset) {
  # Note: pertpop and pertpopconj indicate which population is perturbed by
  # adding bacteria. In case of pertpopconj the added bacteria are
  # plasmid-bearing and replace plasmid-free bacteria. In case of pertpop,
  # the bacteria are plasmid-free and are added to the existing bacteria.
  if(PReplMostAbun == TRUE) {
    pertpop <- "R1"
    pertpopconj <- "P1"
  } else {
    pertpop <- paste0("R", nspecies)
    pertpopconj <- paste0("P", nspecies)
  }
  
  for(abunmodel in abunmodelset) {
    if(abunmodel == "brokenstick") {
      abundance <- brokenstick_fast(nspecies = nspecies, totalabun = totalabun,
                               takelimit = TRUE)
      # Using a number instead of a name, to prevent type-conversion when
      # storing it in matrices.
      abunmodelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt_fast(nspecies = nspecies, totalabun = totalabun,
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
        data <- matrix(data = NA, nrow = nrowdata, ncol = 59 + 4*maxnspecies)
        indexdata <- 1
        relabunRsp <- rep(NA, maxnspecies)
        relabunRconjsp <- rep(NA, maxnspecies)
        relabunPconjsp <- rep(NA, maxnspecies)
        
        for(iter in seq_len(niter)) {
          stableeq <- FALSE
          iterintmat <- 0
          conjratecode <- NA
          taxmatcode <- NA
          
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
          # warning("No stable equilibrium has been found in", niterintmat,
          #         "attempts.")
          # }
          
          # Simulate invasion of plasmid-free bacteria into the plasmid-free
          # equilibrium in the model without conjugation (i.e., test internal
          # stability).
          if(iter == 1L) {
            warning("'tmax' is put at 100 for pertubation with plasmid-free",
                    " bacteria, such that\nequilibrium will be rarely reached.",
                    " This is done to save time and focus on\ninvasion of",
                    " plasmid-bearing bacteria.")
          }
          if(simulateinvasion == TRUE) {
            if(stableeq == FALSE) {
              abunfinal <- perturbequilibrium(abundance = abundance,
                                              intmat = intmat,
                                              growthrate = growthrate,
                                              cost = NULL, conjmat = NULL,
                                              model = "gLV", pertpop = pertpop,
                                              # HARDCODED tmax = 100 to focus on
                                              # invasion of plasmid-bearing
                                              # bacteria, see the warning issued
                                              # for iter == 1L above.
                                              pertmagn = 1, tmax = 100, tstep = 10,
                                              showplot = FALSE, verbose = FALSE,
                                              silentinfgrowth = TRUE,
                                              silenteqnotreached = TRUE)
            } else {
              # No need for simulations if equilibrium is stable
              abunfinal <- list(R = abundance, Rtotal = sum(abundance),
                                npopR = nspecies,
                                P = NULL, Ptotal = NULL, npopP = NULL,
                                pertpopconjsurvived = NULL,
                                timepertpopconjextinct = NULL,
                                timefinal = 1, tmaxshort = 0,
                                eqreached = 1, infgrowth = 0)
            }
          } else {
            # No simulations over time performed, so set values to NA
            abunfinal <- list(R = rep(NA, nspecies), Rtotal = NA,
                              npopR = NA,
                              P = NULL, Ptotal = NULL, npopP = NULL,
                              pertpopconjsurvived = NULL,
                              timepertpopconjextinct = NULL,
                              timefinal = NA, tmaxshort = NA,
                              eqreached = NA, infgrowth = NA)
          }
          
          # Model without conjugation, so Rtotal is total abundance because P
          # does not exist
          relabunRsp[seq_len(nspecies)] <- abunfinal$R / abunfinal$Rtotal
          
          for(cost in costset) {
            conjratecode <- 0
            
            for(conjrate in conjrateset) {
              conjratecode <- conjratecode + 1
              taxmatcode <- 0
              
              for(taxmattype in taxmattypeset) {
                taxmatcode <- taxmatcode + 1
                taxmat <- matrix(rep("SameSpecies", nspecies^2),
                                 nrow = nspecies, ncol = nspecies, byrow = TRUE)
                
                if(taxmattype != "SameSpecies") {
                  # If taxmattype is not 'SameSpecies', 'taxmat' has to be
                  # adjusted to reflect the more-distant relatedness of the
                  # initially plasmid-bearing species to the initially
                  # plasmid-free species. If 'PReplMostAbun' is TRUE the initially
                  # plasmid-bearing species is the most-abundant species (i.e.,
                  # species 1), such that the first row and column of the matrix
                  # taxmat have to be adjusted. If 'PReplMostAbun' is FALSE this
                  # is the least-abundant species (i.e., species nspecies), such
                  # that the last row and column of the matrix taxmat have to be
                  # adjusted. In both cases the diagonal should remain
                  # 'SameSpecies' because those reflect intraspecies
                  # relationships by definition.

                  # Note that PReplMostAbun is used to determine the initially
                  # plasmid-bearing species is also the perturbed
                  # species.
                  if(PReplMostAbun == TRUE) {
                    taxmat[1, -1] <- taxmattype
                    taxmat[-1, 1] <- taxmattype
                  } else {
                    taxmat[nspecies, -nspecies] <- taxmattype
                    taxmat[-nspecies, nspecies] <- taxmattype
                  }
                }
                conjmat <- getconjmat_fast(nspecies = nspecies,
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
                
                # To simulate invasion of plasmid-bearing bacteria into the
                # plasmid-free equilibrium, the abundances of the plasmid-free
                # populations have to be appended to the abundances of the
                # plasmid-bearing populations. 'pertpop' is the population to
                # which bacteria are added.
                if(simulateinvasion == TRUE) {
                  if(eqinfoconj["eigvalRe"] >= 0) {
                    abunfinalconj <- perturbequilibrium(abundance = c(abundance, rep(0, nspecies)),
                                                        intmat = intmat, growthrate = growthrate,
                                                        cost = cost, conjmat = conjmat,
                                                        model = "gLVConj", pertpop = pertpopconj,
                                                        pertmagn = 1, tmax = 5e3, tstep = 10,
                                                        showplot = FALSE, verbose = FALSE,
                                                        silentinfgrowth = TRUE,
                                                        silenteqnotreached = TRUE)
                  } else {
                    # No need for simulations if equilibrium is stable
                    # NOTE: timepertpopconjextinct = 0 is only true if pertmagn is
                    # equal to finalsmallstate. Similarly, timefinal = 1 is not
                    # strictly true, as it will take some time to reach the new
                    # equilibrium
                    abunfinalconj <- list(R = abundance, Rtotal = sum(abundance),
                                          npopR = nspecies,
                                          P = rep(0, nspecies), Ptotal = 0,
                                          npopP = 0, pertpopconjsurvived = 0,
                                          timepertpopconjextinct = 0,
                                          timefinal = 1, tmaxshort = 0,
                                          eqreached = 1, infgrowth = 0)
                  }
                  
                } else {
                  # No simulations over time performed, so set values to NA
                  abunfinalconj <- list(R = rep(NA, nspecies), Rtotal = NA,
                                        npopR = NA,
                                        P = rep(NA, nspecies), Ptotal = NA,
                                        npopP = NA, pertpopconjsurvived = NA,
                                        timepertpopconjextinct = NA,
                                        timefinal = NA, tmaxshort = NA,
                                        eqreached = NA, infgrowth = NA)
                }
                
                # Total (cells / mL) and relative (fractions of total) abundances
                abuntotalconj <- abunfinalconj$Rtotal + abunfinalconj$Ptotal
                relabunRconjsp[seq_len(nspecies)] <- abunfinalconj$R / abuntotalconj
                relabunPconjsp[seq_len(nspecies)] <- abunfinalconj$P / abuntotalconj
                
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
                  cost, conjratecode, taxmatcode, iter, seq_len(nspecies), abundance,
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
                  fracPformedbypertpop, abunfinalconj$pertpopconjsurvived,
                  abunfinalconj$timepertpopconjextinct,
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
                            "eigvalRe", "eigvalIm",
                            "eigvalReSign", "eigvalImSign", "eigvalRep",
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
                            "fracPformedbypertpop", "pertpopconjsurvived",
                            "timepertpopconjextinct",
                            paste0("relabunRsp", seq_len(maxnspecies)),
                            paste0("relabunRconjsp", seq_len(maxnspecies)),
                            paste0("relabunPconjsp", seq_len(maxnspecies)),
                            paste0("relabunconjsp", seq_len(maxnspecies)))
        
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
          # combinations of costs, conjugation rates, and taxmatcodes
          summarydatasimulation <- as_tibble(data) %>%
            group_by(cost, conjratecode, taxmatcode) %>%
            summarise(
              across(c(infgrowth, infgrowthconj, eqreached, eqreachedconj,
                       tmaxshort, tmaxshortconj, pertpopconjsurvived),
                     getfracnotzero, .names = "{.col}{.fn}"),
              timefinalmedian = median(timefinal),
              timefinalconjmedian = median(timefinalconj),
              across(c(abuntotal, abunRtotalconj,
                       starts_with("relabunR"), abunPtotalconj,
                       starts_with("relabunP"), abuntotalconj,
                       relabunRconj, starts_with("relabunconjsp"),
                       starts_with("npop"), nspeciesconj,
                       fracspeciesRconj, fracspeciesPconj, fracPformedbypertpop,
                       timepertpopconjextinct),
                     getsummary4, .names = "{.col}{.fn}"),
              .groups = "drop"
            )
          summarydata <- full_join(summarydata, summarydatasimulation,
                                   by = c("cost", "conjratecode", "taxmatcode"))
        }
        
        rowindexplotdatanew <- rowindexplotdata + length(costset) *
          length(conjrateset) * length(taxmattypeset)
        plotdata[rowindexplotdata:(rowindexplotdatanew - 1), ] <- as.matrix.data.frame(
          tibble(niter, nspecies, abunmodelcode, intmean, selfintmean,
                 summarydata))
        rowindexplotdata <- rowindexplotdatanew
      }
    }
  }
}
duration <- Sys.time() - starttime
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)

colnames(plotdata) <- c("niter", "nspecies", "abunmodelcode", "intmean",
                        "selfintmean",
                        colnames(summarydata))
colnames(datatotal) <- colnames(data)

warnings()
if(simulateinvasion == TRUE) {
  eqnotreached <- 1 - plotdata[, "eqreachedfrac"]
  eqnotreachedconj <- 1 - plotdata[, "eqreachedconjfrac"]
  if(any(eqnotreached > 0)) {
    warning("Fraction of iterations where equilibrium has not been reached",
            " after pertubation with plasmid-free\nbacteria ranges from ",
            signif(summary(eqnotreached)["Min."], 3), " to ",
            signif(summary(eqnotreached)["Max."], 3)," (mean: ",
            signif(summary(eqnotreached)["Mean"], 3), "; median: ",
            signif(summary(eqnotreached)["Median"], 3), "). ",
            " Use silenteqnotreached = FALSE in perturbequilibrium() for more",
            " info.")
  }
  if(any(eqnotreachedconj > 0)) {
    warning("Fraction of iterations where equilibrium has not been reached ",
            " after pertubation with plasmid-bearing\nbacteria ranges from ",
            signif(summary(eqnotreachedconj)["Min."], 3), " to ",
            signif(summary(eqnotreachedconj)["Max."], 3), " (mean: ",
            signif(summary(eqnotreachedconj)["Mean"], 3), "; median: ",
            signif(summary(eqnotreachedconj)["Median"], 3), "). ",
            " Use silenteqnotreached = FALSE in perturbequilibrium() for more",
            " info")
  }
}


#### Saving settings and output to CSV files ####
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
if(PReplMostAbun == FALSE) {
  DateTimeStamp <- paste0(DateTimeStamp, "PReplLeastAbun")
}
if(nrow(plotdata) > 250000) {
  warning("Not saved 'plotdata' to CSV-file because the number of rows (",
          nrow(plotdata), ") exceeds 250000.")
} else {
  write.csv(plotdata, file = paste0(DateTimeStamp, "multispecies.csv"),
            quote = FALSE, row.names = FALSE)
}
names(conjrateset) <- paste0("conjrateset", seq_along(conjrateset))
settings <- c(list(niter = niter, niterintmat = niterintmat,
                   simulateinvasion = simulateinvasion,
                   smallstate = smallstate, smallchange = smallchange,
                   tstep = formals(perturbequilibrium)$tstep,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   costset = costset, conjrateset, taxmattype = taxmattypeset,
                   costtype = costtype,
                   PFrom = if(PReplMostAbun) {"MostAbun"} else {"LeastAbun"},
                   PReplMostAbun = PReplMostAbun, duration = duration))
for(index in seq_along(settings)) {
  # Using write.table instead of write.csv() to be able to use append = TRUE
  write.table(t(as.data.frame(settings[index])), 
              file = paste0(DateTimeStamp, "settings.csv"), append = TRUE,
              quote = FALSE, sep = ",", col.names = FALSE)
}
capture.output(sessionInfo(),
               file = paste0(DateTimeStamp, "sessioninfo_base.txt"))
if(requireNamespace("sessioninfo")) {
  capture.output(sessioninfo::session_info(),
                 file = paste0(DateTimeStamp, "sessioninfo.txt"))
}

if(nrow(datatotal) < 250000) {
  write.csv(datatotal, file = paste0(DateTimeStamp, "multispeciestotal.csv"),
            quote = FALSE, row.names = FALSE)
} else {
  for(nspecies in nspeciesset) {
    filename <- paste0(DateTimeStamp, "multispeciestotal_nsp", nspecies, ".csv")
    write.csv(datatotal[which(near(nspecies, datatotal[, "nspecies"])), ],
              file = filename, quote = FALSE, row.names = FALSE)
    message("Saved file ", filename)
  }
}


#### Reading previously saved data from a CSV file ####
# # To read data from a CSV file, put the CSV file in the working directory
# # (see getwd()), uncomment this section and change the date-time-stamp in the
# # file name (the next line prints a list of files in the working directory
# # that contain 'species' in their name).
# list.files(pattern = "species", ignore.case = TRUE)
# # Note that the extension (.csv) should be included in the file name.
# filename <- file.path("OutputMS", "2023_0428_0507", "2023_0428_0507multispecies.csv")
# plotdata <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                      dec = ".", stringsAsFactors = FALSE)
# # If plotdata has only one column, probably a semicolon instead of a comma was
# # used as separator in the CSV file. So read the file again using that separator.
# if(near(dim(plotdata)[2], 1)) {
#   plotdata <- read.csv(filename, header = TRUE, sep = ";", quote = "\"",
#                        dec = ".", stringsAsFactors = FALSE)
# }
# plotdata <- as.data.frame(plotdata)
# DateTimeStamp <- substr(filename, 1, 16)
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))
# 
# # Similarly for the datatotal CSV file.
# filename <- file.path("OutputMS", "2023_0428_0507",
#                       "2023_0428_0507multispeciestotal_nsp248.csv")
# datatotal <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                      dec = ".", stringsAsFactors = FALSE)
# filename <- file.path("OutputMS", "2023_0428_0507",
#                       "2023_0428_0507multispeciestotal_nsp16.csv")
# datatotal_p2 <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                       dec = ".", stringsAsFactors = FALSE)
# # The first column name of datatotal_p2 is weird, adjust manually
# head(colnames(datatotal))
# head(colnames(datatotal_p2))
# colnames(datatotal_p2)[1] <- colnames(datatotal)[1]
# head(colnames(datatotal_p2))
# datatotal <- rbind(datatotal, datatotal_p2)
# rm(datatotal_p2)


#### Labels and limits for plots ####
labspecies <- paste("Sp.", seq_len(maxnspecies))
names(labspecies) <- seq_len(maxnspecies)
labnspecies <- paste(nspeciesset, "sp.")
names(labnspecies) <- nspeciesset
labmodel <- c("Broken stick", "Dom. preemption")
names(labmodel) <- c(1, 2)
if(bifurparms == FALSE) {
  labcost <- c("Mean fitness\ncost", "High fitness\ncost")
} else {
  labcost <- paste0("Cost: ", costset, "/h")
}
names(labcost) <- costset
labconjrate <- paste("Conjset", seq_along(conjrateset))
names(labconjrate) <- seq_along(conjrateset)
labtaxmat <- c("All conjugation\nrates equal",
               "InitP low inter-\nspecies rates")
names(labtaxmat) <- seq_along(taxmattypeset)
mylabeller <- labeller(species = labspecies, nspecies = labnspecies,
                       abunmodelcode = labmodel,
                       cost = labcost, conjratecode = labconjrate,
                       taxmatcode = labtaxmat, .default = label_value)
plotdata <- as.data.frame(plotdata)
datatotal <- as.data.frame(datatotal)
limitsfraction <- c(0, 1)
# Round the limits to one decimal place, while ensuring that all the data is
# within the rounded limits.
# To do:
# - Switch to safer rounding.
limitsgrowthrate <- c(floor(min(plotdata[, "growthratemin"])*10)/10,
                      ceiling(max(plotdata[, "growthratemax"])*10)/10)
stat_type <- c("min", "mean", "median", "max")
names(stat_type) <- c("Min.", "Mean", "Median", "Max.")


#### To test plots without using CreatePlot() ####
# ggplot(data = plotdata, aes(x = intmean, y = selfintmean, fill = fracstable)) +
#   geom_raster() +
#   theme_bw(base_size = 15) +
#   scale_x_continuous() +
#   scale_y_continuous() +
#   scale_fill_viridis_c("Fraction stable", limits = limitsfraction) +
#   geom_vline(xintercept = 0, col = "grey", size = 1.1) +
#   coord_fixed(ratio = 1, expand = FALSE) +
#   theme(legend.position = "bottom") +
#   labs(x = "Mean interspecies interaction coefficient",
#        y = "Mean intraspecies interaction coefficient",
#        caption = paste(niter, "iterations")) +
#   facet_grid(nspecies ~ abunmodelcode, labeller = mylabeller)


#### Plot output ####
# If the error '$ operator is invalid for atomic vectors' arises, the matrix
# 'plotdata' has not yet been converted to a dataframe: run
# plotdata <- as.data.frame(plotdata) to do so.

if(bifurparms == TRUE) {
  # Add column to dataframe containing conjugation rate
  # NOTE: originally stored as conjratecode because the conjugation rates of the
  # various species can be varied and the variable 'conjrate' is a vector of
  # multiple values, not a single value.
  plotdata$conjrate <- NA
  if(length(seqconjrate) != conjratecode) {
    stop("Value of conjrate code is not equal to length of seqconjrate.
            Conversion of conjratecode to conjrate will be incorrect")
  }
  if(!all(near(conjrate, conjrate[1]))) {
    stop("Species have different conjugation rates, so conversion from",
         " conjugationratecode\nto conjugation rate is not meaningful.\nSee",
         " the settings section in the script and in the file 'settings.csv'",
         " for details")
  }
  for(conjratecode_index in seq_len(conjratecode)) {
    # using dplyr::near() to allow for small (<1e-6) numeric differences
    temp_row_index <- which(near(plotdata[, "conjratecode"], conjratecode_index))
    plotdata[temp_row_index, "conjrate"] <- seqconjrate[conjratecode_index]
  }
  
  
  # ## Show border of ecological stability with heatmap in CreatePlot()
  # # Values in a facet are either all stable, or all unstable (see heatmap above),
  # # making it impossible to plot contours delimiting stable and unstable regions.
  # CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableecol",
  #            filltitle = "fracstableecol", filltype = "continuous",
  #            limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
  #            labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
  #            facetx = "taxmatcode + intmean", facety = "nspecies",
  #            rotate_x_labels = TRUE,
  #            filename = "ecostabxcostyconj")
  # 
  # ## Show border of ecological stability with a contour plot in CreatePlot()
  # # Values in a facet are either all stable, or all unstable (see heatmap above),
  # # making it impossible to plot contours delimiting stable and unstable regions.
  # CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableecol",
  #            contour_var = "fracstableecol", contour_col = "as.factor(nspecies)",
  #            contour_lty = "as.factor(intmean)",
  #            limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
  #            labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
  #            facetx = "taxmatcode", facety = "nspecies",
  #            rotate_x_labels = FALSE, save = FALSE) +
  #   guides(col = guide_legend(ncol = 1), lty = guide_legend(ncol = 1))
  
  
  ## Show border of epidemiological stability with a contour plot in CreatePlot()
  # I set save to FALSE and used ggsave() to ensure the added guides arguments are
  # included in the saved plots.
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
             facetx = "taxmatcode + intmean", facety = "nspecies",
             rotate_x_labels = TRUE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatintmeanynspecies.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(intmean)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
             facetx = "taxmatcode", facety = "nspecies",
             rotate_x_labels = FALSE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspecies.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  # So intmean does not affect the border of stability in conjugation rate/cost-space
  
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
             facetx = "taxmatcode", facety = "intmean",
             rotate_x_labels = TRUE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatyintmean.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  
  # Invasion is possible for very slightly higher costs for a given conjugation
  # rate in a microbiome of 2 species than in a microbiome of 4 or 6 species when
  # the initially plasmid-bearing species belongs to another class.
  
  # When the plasmid-bearing species belongs to another class, invasion is only
  # possible if costs are slightly lower for a given conjugation rate than when
  # all populations belonging to the same species.
  
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "nspecies",
             contour_lty = "as.factor(intmean)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
             facetx = "taxmatcode", facety = ".",
             rotate_x_labels = FALSE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmat.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(taxmatcode)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
             facetx = ".", facety = "nspecies",
             rotate_x_labels = FALSE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabynspecies.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "as.factor(intmean)",
             contour_lty = "as.factor(selfintmean)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(intraspecies conjugation rate of initP)",
             linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
             rotate_x_labels = FALSE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1))
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspeciescolinter.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  
  # Note: assuming sets are chosen such that border of invasion is shown in the
  # plot
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
             contour_var = "fracstableepi", contour_col = "as.factor(nspecies)",
             contour_lty = "as.factor(taxmatcode)",
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
             facetx = ".", facety = ".",
             rotate_x_labels = FALSE, save = FALSE) +
    theme(legend.box = "vertical", legend.margin = margin(rep(-5, 4), unit = "pt")) +
    guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
    geom_point(data = expand.grid(cost = c(0.05, 0.09), conjrate = -12),
               inherit.aes = FALSE, mapping = aes(x = cost, y = conjrate),
               color = "black", size = 2) +
    annotate("text", label = c("Invasion of plasmid-\nbearing bacteria",
                               "No invasion of plasmid-\nbearing bacteria"),
             x = quantile(costset, c(0.1, 0.9)),
             y = quantile(log10(seqconjrate), c(0.9, 0.1)), 
             hjust = "inward", vjust = "inward", size = 4)
  if(saveplots == TRUE) {
    ggsave(paste0(DateTimeStamp, "epistab.png"),
           width = 2150, height = 2150, units = "px", dpi = 300)
  }
  
  # Need to set filltype to continuous to prevent error on missing filllabels
  CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableepi",
             filltitle = "fracstableepi", contour_var = NULL, contour_col = NULL,
             contour_lty = NULL, filltype = "continuous", limx = range(costset),
             limx = c(0, NA), limy = range(log10(seqconjrate)), ratio = NULL,
             title = "Epidemiological (in)stability",
             labx = "Cost", laby = "Log10(conjugation rate)", linezero = FALSE,
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
CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "fracstableecol",
           filltitle = "Fraction ecologically\nstable",
           filltype = "continuous", limits = limitsfraction, tag = "A",
           filename = "Fig09A", width = 1650, height = 2675)
CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "fracstableecol",
           filltitle = "Fraction ecologically\nstable",
           filltype = "continuous", limits = limitsfraction, tag = "A",
           filename = "Fig09A_altres", width = 825, height = 2675)
CreatePlot(fillvar = "fracstableepi",
           filltitle = "Fraction epidemiologically\nstable",
           filltype = "continuous", limits = limitsfraction)

# Show the effect of adding conjugation on stability
CreatePlot(fillvar = "fracunstableunstable + fracneutralneutral + fracstablestable",
           filltitle = "Fraction stability\nunchanged\nby conjugation",
           filltype = "continuous", limits = limitsfraction)

eq_states <- c("unstable", "neutral", "stable")
for(eq_status_without in eq_states) {
  for(eq_status_with in eq_states) {
    print(CreatePlot(fillvar = paste0("frac", eq_status_without, eq_status_with),
                     filltitle = paste0("Fraction ",
                                        eq_status_without, " without and\n",
                                        eq_status_with, " with conjugation"),
                     filltype = "continuous", limits = limitsfraction))
  }
}

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
  ggsave(paste0(DateTimeStamp, "eigenvaluesdistr.png"),
         width = 1650, height = 2675, units = "px", dpi = 300)
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
  ggsave(paste0(DateTimeStamp, "eigenvaluesdistrconj.png"),
         width = 1650, height = 2675, units = "px", dpi = 300)
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
  CreatePlot(fillvar = "infgrowthconjfrac",
             filltitle = "Fraction infinite growth\nwith conjugation",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             filename = "Fig13A")
  ggplot(data = plotdata, aes(x = intmean, y = selfintmean,
                              fill = infgrowthconjfrac < 0.75 | intmean > 5e-12)) +
    geom_raster() +
    theme_bw(base_size = 11) +
    scale_x_continuous() +
    scale_y_continuous() +
    scale_fill_viridis_d(limits = NULL) +
    coord_fixed(ratio = 1, expand = FALSE) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "Mean interspecies interaction coefficient",
         y = "Mean intraspecies interaction coefficient",
         caption = paste(niter, "iterations"),
         tag = "A") +
    facet_grid(conjratecode + nspecies ~ cost + taxmatcode,
               labeller = label_both, as.table = TRUE, drop = TRUE)
  if(saveplots == TRUE) {
    ggsave("Fig10A_alt.png")
  }
  
  CreatePlot(fillvar = "eqreachedconjfrac",
             filltitle = "Fraction equilibrium\nreached with\nconjugation",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedconjfrac",
             filltitle = "Fraction equilibrium\nreached with\nconjugation",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             filename = "Fig12A")
  CreatePlot(fillvar = "tmaxshortconjfrac",
             filltitle = "Fraction tmax too short\nwith conjugation",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "tmaxshortconjfrac",
             filltitle = "Fraction tmax too short\nwith conjugation",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             filename = "Fig14A")
}
CreatePlot(fillvar = "fracreal", filltitle = "Fraction real\neigenvalues",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracrealconj",
           filltitle = "Fraction real\neigenvalues\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fraceigvalRep", filltitle = "Fraction repeated\neigenvalues",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fraceigvalRepconj",
           filltitle = "Fraction repeated\neigenvalues\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)

## Growth rates
# Plot summary data for the calculated growth rates 
for(ind_stat_type in seq_along(stat_type)) {
  print(CreatePlot(fillvar = paste0("growthrate", stat_type[ind_stat_type]),
                   filltitle = paste(names(stat_type[ind_stat_type]), "growth rate"),
                   filltype = "continuous", limits = limitsgrowthrate))
}

# Plot mean and median growth rate with their own limits
CreatePlot(fillvar = "growthratemean",
           filltitle = "Mean growth rate",
           filltype = "continuous", limits = NULL,
           filename = "multispeciesgrowthratemeancontinuous_ownlimits")
CreatePlot(fillvar = "growthratemedian",
           filltitle = "Median growth rate",
           filltype = "continuous", limits = NULL,
           filename = "multispeciesgrowthratemediancontinuous_ownlimits")

# Show the relation of interactions and species-specific growth rate required to
# obtain an equilibrium. Costs and conjugation rate do not affect growth rate,
# so data has been filtered to have only one value for them.
datatotalfiltercostconj <- filter(datatotal, near(cost, costset[1]),
                                  near(conjratecode, 1), near(taxmatcode, 1))
datatotalfiltercostconjsp1 <- filter(datatotalfiltercostconj,
                                     near(species, 1))
ggplot(data = datatotalfiltercostconjsp1,
  aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_point(aes(col = iter), size = 0.1) +
  facet_grid(rows = vars(selfintmean), cols = vars(nspecies),
             labeller = mylabeller, drop = TRUE, scales = "fixed") +
  scale_color_viridis_c() +
  labs(title = "Growth rate of species 1",
       caption = paste(niter, "iterations")) +
  guides(color = guide_colourbar())
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsintmean_sp1_col.png"),
         width = 3300, height = 2675, units = "px", dpi = 300)
}

ggplot(data = datatotalfiltercostconjsp1,
       aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_point(size = 0.1) +
  facet_grid(rows = vars(selfintmean), cols = vars(nspecies),
             labeller = mylabeller, drop = TRUE, scales = "fixed") +
  scale_color_viridis_c() +
  labs(title = "Growth rate of species 1",
       caption = paste(niter, "iterations")) +
  guides(color = guide_colourbar())
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsintmean_sp1_black.png"),
         width = 3300, height = 2675, units = "px", dpi = 300)
}

datatotalfiltercostconj_iter1 <- datatotalfiltercostconj[
  which(near(1, datatotalfiltercostconj[, "iter"])), ]

# To do:
# - The colorbar labels should be shifted to the left
ggplot(data = datatotalfiltercostconj, aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = selfintmean), size = 1) +
  facet_grid(rows = vars(species, abunmodelcode), cols = vars(nspecies),
             labeller = mylabeller, drop = TRUE) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations")) +
  guides(color = guide_colourbar(label.hjust = 1, label.vjust = 1,
                                 label.theme = element_text(angle = 45)))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsintmean.png"),
         width = 3300, height = 2675, units = "px", dpi = 300)
}

ggplot(data = datatotalfiltercostconj, aes(x = selfintmean, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = intmean), size = 1) +
  facet_grid(rows = vars(nspecies), cols = vars(species, abunmodelcode),
             labeller = mylabeller) +
  scale_color_viridis_c() +
  labs(caption = paste(niter, "iterations")) +
  guides(color = guide_colourbar(label.hjust = 1, label.vjust = 1,
                                 label.theme = element_text(angle = 45)))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsselfintmean.png"),
         width = 3300, height = 2675, units = "px", dpi = 300)
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
for(ind_stat_type in seq_along(stat_type)) {
  print(CreatePlot(fillvar = paste0("iterintmat", stat_type[ind_stat_type]),
                   filltitle = paste(names(stat_type[ind_stat_type]),
                                     "number of\niterations to reach\nstable equilibrium"),
                   filltype = "continuous", limits = c(1, niterintmat)))
}

if(simulateinvasion == TRUE) {
  subplasmidfree <- "Perturbation with plasmid-free bacteria"
  subplasmidbearing <- "Perturbation with plasmid-bearing bacteria"
  if(PReplMostAbun == FALSE) {
    title_add <- " (PReplLeastAbun)"
  } else {
    title_add <- NULL
  }
  
  limitstime <- range(plotdata[, "timefinalmedian"],
                      plotdata[, "timefinalconjmedian"])
  title <- paste0("Time to reach equilibrium after perturbation", title_add)
  CreatePlot(fillvar = "timefinalmedian", filltitle = "Median time",
             filltype = "continuous", limits = limitstime,
             title = title, subtitle = subplasmidfree, rotate_legend = TRUE)
  CreatePlot(fillvar = "timefinalconjmedian", filltitle = "Median time",
             filltype = "continuous", limits = limitstime,
             title = title, subtitle = subplasmidbearing, rotate_legend = TRUE)
  
  ## Plots on survival and extinction after perturbation
  limitsmeannspecies <- range(plotdata[, "npopRmean"],
                              plotdata[, "nspeciesconjmean"], finite = TRUE)
  title <- paste0("Number of species surviving after perturbation", title_add)
  # Note: In the model without plasmids, the number of populations is equal to
  # the number of species.
  CreatePlot(fillvar = "npopRmean", filltitle = "Mean total number\nof species",
             filltype = "continuous", limits = limitsmeannspecies,
             title = title, subtitle = subplasmidfree,
             filename = "nspeciesRmean")
  CreatePlot(fillvar = "nspeciesconjmean", filltitle = "Mean total number\nof species",
             filltype = "continuous", limits = limitsmeannspecies,
             title = title, subtitle = subplasmidbearing)
  
  title <- paste0("Number of species extinct after perturbation", title_add)
  CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                               near(conjratecode, 1), near(taxmatcode, 1)),
             fillvar = "nspecies - npopRmean",
             filltitle = "Mean number of\nspecies extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidfree,
             facetx = "abunmodelcode", facety = "nspecies",
             filename = "nspeciesRmeanextinctfewfacets")
  CreatePlot(fillvar = "nspecies - npopRmean",
             filltitle = "Mean number of\nspecies extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidfree,
             filename = "nspeciesRmeanextinct")
  CreatePlot(fillvar = "nspecies - nspeciesconjmean",
             filltitle = "Mean number of\nspecies extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidbearing,
             filename = "nspeciesconjmeanextinct")
  
  title <- paste0("Fraction of species extinct after perturbation", title_add)
  CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                               near(conjratecode, 1), near(taxmatcode, 1)),
             fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of\nspecies extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree,
             facetx = "abunmodelcode", facety = "nspecies",
             filename = "fracspeciesextinctmeanfewfacets")
  CreatePlot(fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of\nspecies extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree,
             filename = "fracspeciesextinctmean")
  CreatePlot(fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of\nspecies extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree)
  CreatePlot(fillvar = "1 - (nspeciesconjmean / nspecies)",
             filltitle = "Mean fraction of\nspecies extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing,
             filename = "fracspeciesextinctmeanconj")
  CreatePlot(fillvar = "1 - (nspeciesconjmean / nspecies)",
             filltitle = "Mean fraction of\nspecies extinct",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             title = title, subtitle = subplasmidbearing,
             filename = "Fig10A")
  
  title <- paste0("Fraction of species after perturbation", title_add)
  CreatePlot(fillvar = "fracspeciesRconjmean",
             filltitle = "Mean fraction of surviving species\nwith a plasmid-free population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "npopRconjmean / nspecies",
             filltitle = "Mean fraction of initial species\nwith a plasmid-free population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "fracspeciesPconjmean",
             filltitle = "Mean fraction of surviving species\nwith a plasmid-bearing population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "fracspeciesPconjmean",
             filltitle = "Mean fraction of surviving species\nwith a plasmid-bearing population",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             title = title, subtitle = subplasmidbearing,
             filename = "Fig07A")
  CreatePlot(fillvar = "npopPconjmean / nspecies",
             filltitle = "Mean fraction of initial species\nwith a plasmid-bearing population",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "npopPconjmean / nspecies",
             filltitle = "Mean fraction of initial species\nwith a plasmid-bearing population",
             filltype = "continuous", limits = limitsfraction, tag = "B",
             title = title, subtitle = subplasmidbearing, filename = "Fig10B")
  
  ## Plots of fractions of bacteria that are plasmid-free or plasmid-bearing
  # after perturbations. Only abundances where equilibrium was reached are
  # considered.
  title <- paste0("Fraction bacteria after perturbation", title_add)
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("relabunRconj", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of bacteria\nthat is plasmid-free"),
                     filltype = "continuous", limits = limitsfraction,
                     title = title, subtitle = subplasmidbearing))
    
    print(CreatePlot(fillvar =  paste0("relabunPconj", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of bacteria\nthat is plasmid-bearing"),
                     filltype = "continuous", limits = limitsfraction,
                     title = title, subtitle = subplasmidbearing))
    
    print(CreatePlot(fillvar =  paste0("fracPformedbypertpop", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of plasmid-bearing\nbacteria",
                                       "belonging to the\ninitially plasmid-bearing species"),
                     filltype = "continuous", limits = limitsfraction))
    
    print(CreatePlot(fillvar =  paste0("timepertpopconjextinct", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "time initially plasmid-\nbearing",
                                       "population went\nextinct"),
                     filltype = "continuous", rotate_legend = TRUE))
  }
  
  CreatePlot(fillvar =  "relabunPconjmean",
             filltitle = "Mean fraction of bacteria\nthat is plasmid-bearing",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             title = title, subtitle = subplasmidbearing, filename = "Fig03A")
  
  CreatePlot(fillvar =  "1 - fracPformedbypertpopmean",
             filltitle = paste("Mean fraction of plasmid-bearing\nbacteria",
                               "belonging to the\ninitially plasmid-free",
                               "species"),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             filename = "Fig04B")
  
  CreatePlot(fillvar =  "fracPformedbypertpopmean",
             filltitle = paste("Mean fraction of plasmid-bearing\nbacteria",
                               "belonging to the\ninitially plasmid-bearing species"),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             filename = "Fig04B_alt")
  
  CreatePlot(fillvar = "pertpopconjsurvivedfrac",
             filltitle = paste("Fraction of iterations where\ninitially",
                               "plasmid-bearing\npopulation survived"),
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "pertpopconjsurvivedfrac",
             filltitle = paste("Fraction of iterations where\ninitially",
                               "plasmid-bearing\npopulation survived"),
             filltype = "continuous", limits = limitsfraction, tag = "A",
             filename = "Fig06A")
  
  ## Plot total abundances after perturbation with plasmid-free bacteria in
  # models without plasmids. Only abundances where equilibrium was reached are
  # considered. Although costs and conjugation rates do not influence the
  # outcome, they are included as facets to facilitate comparison with plots of
  # abundances after perturbation with plasmid-bearing bacteria.
  
  # log10(0) leads to -Inf and is displayed as a gray square in the plots,
  # leading to different numbers of gray squares between plots. It is also
  # confusing because gray squares can also indicate no data is available
  # because all simulations were discarded because equilibrium was not reached.
  # Therefore I plot log10(1 + x) instead of log10(x). I do not use the
  # built-in function log1p(x) because that uses natural logarithms. For
  # relative abundances I used log(1e-6 + x) because adding 1 is too much.
  title <- paste0("Total abundance after perturbation", title_add)
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abuntotal", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "abundance of plasmid-\nfree bacteria)"),
                     filltype = "continuous", title = title, subtitle = subplasmidfree))
  }
  
  ## Plot of total abundances after perturbation with plasmid-bearing bacteria in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abuntotalconj", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "total\nabundance)"),
                     filltype = "continuous", title = title, subtitle = subplasmidbearing))
  }
  
  ## Plot total abundances of plasmid-free populations after perturbations in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abunRtotalconj", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "abundance of plasmid-\nfree bacteria)"),
                     filltype = "continuous", title = title, subtitle = subplasmidbearing))
  }
  
  ## Plot total abundances of plasmid-bearing populations after perturbations for
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abunPtotalconj", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "abundance of plasmid-\nbearing bacteria)"),
                     filltype = "continuous", title = title, subtitle = subplasmidbearing))
  }
  
  ## Plots comparing species abundances after perturbation with plasmid-bearing
  # bacteria
  if(PReplMostAbun == TRUE) {
    add_filltitle <- "after\nadding some R of the most-abundant sp."
    add_filltitleconj <- "after\nreplacing some R of the most-abundant\nsp. with P of that sp."
  } else {
    add_filltitle <- "after\nadding some R of the least-abundant sp."
    add_filltitleconj <- "after\nreplacing some R of the least-abundant\nsp. with P of that sp."
  }
  
  limits <- range(c(plotdata[, "relabunRsp1mean"],
                    plotdata[, "relabunconjsp1mean"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp1mean",
             filltitle = paste("Mean rel. abundance of sp1", add_filltitle),
             filltype = "continuous", limits = limits, rotate_legend = TRUE,
             filename = "relabunRsp1meancontinuouschangedlim")
  CreatePlot(fillvar = "relabunconjsp1mean",
             filltitle = paste("Mean rel. abundance of sp1", add_filltitleconj),
             filltype = "continuous", limits = limits, rotate_legend = TRUE,
             filename = "relabunconjsp1meancontinuouschangedlim")
  
  limits <- range(c(plotdata[, "relabunRsp1median"],
                    plotdata[, "relabunconjsp1median"]), na.rm = TRUE)
  CreatePlot(fillvar = "relabunRsp1median",
             filltitle = paste("Median rel. abundance of sp1", add_filltitle),
             filltype = "continuous", limits = limits, rotate_legend = TRUE,
             filename = "relabunRsp1mediancontinuouschangedlim")
  CreatePlot(fillvar = "relabunconjsp1median",
             filltitle = paste("Median rel. abundance of sp1", add_filltitleconj),
             filltype = "continuous", limits = limits, rotate_legend = TRUE,
             filename = "relabunconjsp1mediancontinuouschangedlim")
  
  # Print relative abundance of the first and last species after perturbation
  # without and with plasmids, on a normal scale and a log scale.
  for(species_i in unique(c(1, nspeciesset)))  {
    for(ind_stat_type in seq_along(stat_type)) {
      print(CreatePlot(fillvar = paste0("relabunRsp", species_i,
                                        stat_type[ind_stat_type]),
                       filltitle = paste0(names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitle),
                       filltype = "continuous", limits = limitsfraction))
      
      print(CreatePlot(fillvar = paste0("log10(1e-6 + relabunRsp", species_i,
                                        stat_type[ind_stat_type], ")"),
                       filltitle = paste0("Log10(1e-6 + ",
                                          names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitle, ")"),
                       filltype = "continuous", rotate_legend = TRUE))
      
      print(CreatePlot(fillvar = paste0("relabunconjsp", species_i,
                                        stat_type[ind_stat_type]),
                       filltitle = paste0(names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitleconj),
                       filltype = "continuous", limits = limitsfraction))
      
      print(CreatePlot(fillvar = paste0("log10(1e-6 + relabunconjsp", species_i,
                                        stat_type[ind_stat_type], ")"),
                       filltitle = paste0("Log10(1e-6 + ",
                                          names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitleconj, ")"),
                       filltype = "continuous", rotate_legend = TRUE))
    }
  }
}

CreatePlot(fillvar = "relabunconjsp1mean",
           filltitle = paste0("Mean rel. abundance of initially",
                              "\nplasmid-bearing species ", add_filltitleconj),
           filltype = "continuous", limits = limitsfraction, tag = "A",
           filename = "Fig04A")

CreatePlot(fillvar = "log10(1e-6 + relabunconjsp1mean)",
           filltitle = paste0("log10(1e-6 + Mean rel. abundance of initially",
                              "\nplasmid-bearing species) ", add_filltitleconj),
           filltype = "continuous", limits = NULL, tag = "A",
           filename = "Fig04A_log")

CreatePlot(fillvar = "1 - relabunconjsp1mean",
           filltitle = paste0("Mean rel. abundance of initially",
                              "\nplasmid-free species ", add_filltitleconj),
           filltype = "continuous", limits = limitsfraction, tag = "A",
           filename = "Fig04A_alt")

CreatePlot(fillvar = "log10(1 - relabunconjsp1mean + 1e-6)",
           filltitle = paste0("log10(1e-6 + Mean rel. abundance of initially",
                              "\nplasmid-free species) ", add_filltitleconj),
           filltype = "continuous", limits = NULL, tag = "A",
           filename = "Fig04A_alt_log")

#### Compare abundance models ####
compareabun <- NULL
for(abunmodel in abunmodelset) {
  for(nspecies in nspeciesset) {
    comparetemp <- data.frame(
      nspecies = as.factor(nspecies),
      species = seq_len(nspecies),
      abun = switch(abunmodel,
                    "brokenstick" = brokenstick(nspecies = nspecies,
                                                totalabun = totalabun,
                                                takelimit = TRUE),
                    "dompreempt" = dompreempt(nspecies = nspecies,
                                              totalabun = totalabun,
                                              takelimit = TRUE),
                    stop("'abunmodel' should be 'brokenstick' or 'dompreempt'.")),
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
  scale_x_continuous(limits = c(1, maxnspecies)) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_grid(rows = facet_rows, cols = vars(nspecies), labeller = mylabeller) +
  theme(panel.grid = element_line(size = 1),
        legend.position = "bottom", legend.justification = "left",
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
  scale_x_continuous(limits = c(1, maxnspecies)) +
  scale_y_continuous(trans = "log10") +
  facet_grid(rows = facet_rows, cols = vars(nspecies), labeller = mylabeller) +
  theme(panel.grid = element_line(size = 1),
        legend.position = "bottom", legend.justification = "left",
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
ggsave(filename = paste0(DateTimeStamp, "Fig08.png"), plot = plotcompareabunlog,
       width = 16, height = 16, units = "cm")

write.csv(compareabun,
          file = paste0(DateTimeStamp, "compareabunmodel.csv"),
          quote = FALSE, row.names = FALSE)

for(nspecies in nspeciesset) {
  abundata_temp <- compareabun[which(compareabun[, "nspecies"] %in% nspecies),
                               c("species", "abun")]
  print(
    ggplot(data = abundata_temp, aes(x = species, y = abun)) +
      theme_bw(base_size = 14) +
      scale_x_continuous(limits = c(1, NA), expand = c(0.025, 0.025)) +
      scale_y_continuous(limits = c(0, 8e10), expand = c(0.025, 0.025)) +
      theme(panel.grid = element_line(size = 1),
            legend.position = "bottom", legend.justification = "left",
            legend.margin = margin(-5, 0, -5, 0),
            legend.box = "vertical", legend.box.just = "left",
            legend.box.margin = margin(-5, 0, 0, 0)) +
      labs(title = "Example species abundances", subtitle = "Linear scale",
           x = "Species rank", y = "Species abundance") +
      geom_line(size = 1.25) +
      geom_point(size = 3)
  )
  if(saveplots == TRUE) {
    filename <- paste0(DateTimeStamp, "exampleabun_", nspecies, "sp.png")
    ggsave(filename, width = 16, height = 16, units = "cm")
  }
  
  print(
    ggplot(data = abundata_temp, aes(x = species, y = abun)) +
      theme_bw(base_size = 14) +
      scale_x_continuous(limits = c(1, NA), expand = c(0.025, 0.025)) +
      scale_y_continuous(trans = "log10") +
      theme(panel.grid = element_line(size = 1),
            legend.position = "bottom", legend.justification = "left",
            legend.margin = margin(-5, 0, -5, 0),
            legend.box = "vertical", legend.box.just = "left",
            legend.box.margin = margin(-5, 0, 0, 0)) +
      labs(title = "Example species abundances", subtitle = "Logarithmic scale",
           x = "Species rank", y = "Species abundance") +
      geom_line(size = 1.25) +
      geom_point(size = 3)
  )
  if(saveplots == TRUE) {
    filename <- paste0(DateTimeStamp, "exampleabunlog_", nspecies, "sp.png")
    ggsave(filename, width = 16, height = 16, units = "cm")
  }
}
