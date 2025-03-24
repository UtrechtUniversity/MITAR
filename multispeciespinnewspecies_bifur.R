################################################################################
## Modelling the effects of ecological interactions and distinct conjugation  ##
## rates on the invasion of a conjugative plasmid in bacterial communities    ##
################################################################################


#### Introduction ####
# Use bifurcation-like plots of the result of simulations with a generalised
# Lotka-Volterra model elaborated with plasmid-bearing populations and
# conjugation to show how the intraspecies conjugation rate of the initially
# plasmid-bearing species should change to enable invasion of plasmids if the
# interspecies conjugation rate to and from that species are reduced. The
# bifurcation-like plots in the main scripts reflect a different scenario: they
# show how all conjugation rates simultaneously should change to enable invasion
# of plasmids. 

# NOTE:
# - Some functions differ from the functions in the script multispecies_bifur.R
#   which models the case where the plasmid is introduced through an
#   already-present species.


#### References ####
# Alderliesten JB, Duxbury SJN, Zwart MP, de Visser JAGM, Stegeman S, Fischer
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


#### Loading required packages ####
library(dplyr)     # across(), full_join(), group_by(), near(), summarise()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()
# On the pipe operator (%>%), see ?'%>%' and Ch. 18 'Pipes' in Wickham 2017.


#### Settings and defining parameter space ####
# Note: to simulate that each species belongs to a different class, an
# additional taxmattype has to be added. Then only the interspecies conjugation
# rates should be reduced, as the intraspecies conjugation rates of all species
# are still those given in conjrateset. This unchanged intraspecies conjugation
# rate makes it different from reducing conjrateset 1000-fold and using
# 'taxmatsame' as taxmattype.

## Parameter set to create bifurcation-like plots showing the border of
# epidemiological stability in the conjugation rate/cost space
saveplots <- TRUE
niterintmat <- 1
smallstate <- NA
finalsmallstate <- NA
smallchange <- NA
totalabun <- 1e11
nspeciesset <- 1 + c(2, 4, 8, 16) # Because species do matter for threshold of costs
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
newgrowthratecode <- 2
costset <- seq(from = 0, to = 1, by = 0.0025)
costtype <- "absolute"
costmark <- NULL # Plot dotted vertical lines at indicated values if not NULL
conjrate_base <- 1e-12
seqconjrate <- 10^seq(from = -13.0, to = -11.5, by = 0.05)
# If taxmattype is "SameSpecies", the conjugation rate is the same for all
# populations, and equal to 'conjrateset' given above. If taxmattype is
# "OtherClass", the interspecies conjugation rate to and from the initially
# plasmid-bearing population (the newly added species 1) is reduced a 1000-fold.
taxmattypeset <- c("SameSpecies", "OtherClass")
# To plot 16 species need 16 colours, currently only 11 so repeat them.
mycol <- rep(c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
               "darkorange", "green1", "yellow", "hotpink"), 2)

# The conjugation rate given here is the 'overall' conjugation rate. For all
# species, the intraspecies conjugation rate will be equal to this overall
# conjugation rate. The interspecies conjugation rate will be equal to this
# overall conjugation rate if all populations belong to the same species.
# If the initially plasmid-bearing species is considered to be from another
# taxonomic class than the initially plasmid-free species, these overall
# conjugation rates are reduced a thousand-fold to obtain the interspecies
# conjugation rate between initially plasmid-bearing and initially plasmid-free
# species.
# See also 'conjrate_base' defined above.
conjrateset <- NULL
for(conjrate in seqconjrate) {
  conjrateset <- c(conjrateset, list(rep(conjrate, maxnspecies)))
}
niter <- 1
intmeanset <- c(-1e-11, 0, 5e-12)
selfintmeanset <- c(-1e-11, -5e-12, 0)


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
                       marginx = NULL, marginy = NULL, base_size = 13,
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
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = 3*4 + 30)
nrowdatatotal <- prod(lengths(list(abunmodelset,intmeanset, selfintmeanset,
                                   costset, conjrateset, taxmattypeset),
                              use.names = FALSE))*niter*sum(nspeciesset)
indexdatatotal <- 1

# Run simulations
rowindexplotdata <- 1
rowindexdata <- 1

for(nspecies in nspeciesset) {
  for(abunmodel in abunmodelset) {
    if(abunmodel == "brokenstick") {
      abundance <- brokenstick_fast(nspecies = nspecies - 1, totalabun = totalabun,
                                    takelimit = TRUE)
      # Using a number instead of a name, to prevent type-conversion when
      # storing it in matrices.
      abunmodelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt_fast(nspecies = nspecies - 1, totalabun = totalabun,
                                   takelimit = TRUE)
      abunmodelcode <- 2
    }
    
    for(intmean in intmeanset) {
      
      for(selfintmean in selfintmeanset) {
        print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
                     ", intmean = ", intmean, ", selfintmean = ", selfintmean,
                     ": started at ", Sys.time()), quote = FALSE)
        nrowdata <- niter * nspecies * length(costset) * length(conjrateset) *
          length(taxmattypeset)
        data <- matrix(data = NA, nrow = nrowdata, ncol = 36)
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
            growthrateeq <- getgrowthrate(abundance = abundance,
                                          intmat = intmat[-1, -1])
            # The growth rate of new species is the mean growth rate of the
            # plasmid-free species decreased with 2 standard deviations,
            # unchanged, or increased with 2 standard deviations when
            # 'newgrowthratecode' is 1, 2, or 3, respectively.
            growthrate <- c(switch(newgrowthratecode,
                                   mean(growthrateeq) - 2*sd(growthrateeq),
                                   mean(growthrateeq),
                                   mean(growthrateeq) + 2*sd(growthrateeq)),
                            growthrateeq)
            eqinfo <- geteqinfo(model = "gLV", abundance = c(0, abundance),
                                intmat = intmat, growthrate = growthrate)
            if(eqinfo["eigvalRe"] < 0) {
              stableeq <- TRUE
            }
            iterintmat <- iterintmat + 1
          }
          # To get a warning if the plasmid-free equilibrium is not stable,
          # uncomment next lines.
          # if(stableeq == FALSE) {
          #   warning("No stable equilibrium has been found in", niterintmat,
          #           " attempts.")
          # }
          
          for(cost in costset) {
            if(abs(cost %% (costset[2] * 5)) < 1e-5) {
              print(paste0("cost = ", cost, ", started at ", Sys.time()),
                    quote = FALSE)
            }
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
                  # initially plasmid-bearing species to the initially plasmid-
                  # free species. In this model the initially plasmid-bearing
                  # species is always the newly-introduced species 1, such that
                  # the first row and column of the matrix taxmat have to be
                  # adjusted. The diagonal should remain 'SameSpecies' because
                  # those reflect intraspecies relationships by definition.
                  taxmat[1, -1] <- taxmattype
                  taxmat[-1, 1] <- taxmattype
                }
                conjmat <- getconjmat_fast(nspecies = nspecies,
                                           conjrate = conjrate, taxmat = taxmat)
                conjmat[-1, -1] <- conjrate_base
                if(taxmattype == "SameSpecies") {
                  conjmat[1, ] <- conjrate_base
                  conjmat[, 1] <- conjrate_base
                } else {
                  conjmat[1, ] <- conjrate_base / 1000
                  conjmat[, 1] <- conjrate_base / 1000
                }
                conjmat[1, 1] <- conjrate[1]
                
                # Get equilibrium characteristics for plasmid-free equilibrium in
                # the model with conjugation
                eqinfoconj <- geteqinfo(model = "gLVConj",
                                        abundance = c(0, abundance, rep(0, nspecies)),
                                        intmat = intmat, growthrate = growthrate,
                                        cost = cost, conjmat = conjmat)
                
                # Get equilibrium characteristics regarding ecological and
                # epidemiological stability of the plasmid-free equilibrium (see
                # Roberts and Heesterbeek 2021). The abundance of 0 is the
                # initially plasmid-bearing population.
                eqinfoecol <- geteqinfo(model = "ecol",
                                        abundance = c(0, abundance, rep(0, nspecies)),
                                        intmat = intmat, growthrate = growthrate,
                                        cost = cost, conjmat = conjmat)
                eqinfoepi <- geteqinfo(model = "epi",
                                       abundance = c(0, abundance, rep(0, nspecies)),
                                       intmat = intmat, growthrate = growthrate,
                                       cost = cost, conjmat = conjmat)
                
                # Append the signs of the real parts of the largest eigenvalues to
                # indicate change in stability without and with the plasmid (-1 =
                # stable, 0 = neutral, 1 = unstable). The + 3 ensure numbers are
                # positive, to prevent paste0(-1, -1) leading to NA because of
                # as.integer(-1-1))
                compstability <- as.integer(paste0(eqinfo["eigvalReSign"] + 3,
                                                   eqinfoconj["eigvalReSign"] + 3))
                
                nspeciesconj <- NA 
                indexdatanew <- indexdata + nspecies
                
                data[indexdata:(indexdatanew - 1), ] <- cbind(
                  niter, nspecies, abunmodelcode, intmean, selfintmean,
                  cost, conjratecode, taxmatcode, iter, seq_len(nspecies), c(0, abundance),
                  diag(intmat), c(growthrate), newgrowthratecode, iterintmat,
                  matrix(rep(eqinfo, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(eqinfoconj, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(eqinfoecol, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(eqinfoepi, nspecies), nrow = nspecies, byrow = TRUE),
                  compstability)
                indexdata <- indexdatanew
              }
            }
          }
        }
        indexdatatotal <- indexdatatotal + nrowdata
        
        colnames(data) <- c("niter", "nspecies", "abunmodelcode",
                            "intmean", "selfintmean", "cost", "conjratecode",
                            "taxmatcode", "iter", "species", "abundance",
                            "selfintdata", "growthrate", "newgrowthratecode",
                            "iterintmat",
                            "eigvalRe", "eigvalIm",
                            "eigvalReSign", "eigvalImSign", "eigvalRep",
                            "eigvalReconj", "eigvalImconj",
                            "eigvalReSignconj", "eigvalImSignconj", "eigvalRepconj",
                            "eigvalReecol", "eigvalImecol",
                            "eigvalReSignecol", "eigvalImSignecol", "eigvalRepecol",
                            "eigvalReepi", "eigvalImepi",
                            "eigvalReSignepi", "eigvalImSignepi", "eigvalRepepi",
                            "compstability")
        
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
        
        rowindexplotdatanew <- rowindexplotdata + length(costset) *
          length(conjrateset) * length(taxmattypeset)
        plotdata[rowindexplotdata:(rowindexplotdatanew - 1), ] <- as.matrix.data.frame(
          tibble(niter, nspecies, abunmodelcode, intmean, selfintmean,
                 newgrowthratecode, summarydata))
        rowindexplotdata <- rowindexplotdatanew
      }
    }
  }
}
duration <- Sys.time() - starttime
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)

colnames(plotdata) <- c("niter", "nspecies", "abunmodelcode", "intmean",
                        "selfintmean", "newgrowthratecode",
                        colnames(summarydata))
summary(warnings())
rm(summarydata)


#### Saving settings and output to CSV files ####
DateTimeStamp <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M"), "PInNewSp")
if(nrow(plotdata) > 250000) {
  warning("Not saved 'plotdata' to CSV-file because the number of rows (",
          nrow(plotdata), ") exceeds 250000.")
} else {
  write.csv(plotdata, file = paste0(DateTimeStamp, "multispecies.csv"),
            quote = FALSE, row.names = FALSE)
}
names(conjrateset) <- paste0("conjrateset", seq_along(conjrateset))
settings <- c(list(niter = niter, niterintmat = niterintmat,
                   simulateinvasion = FALSE,
                   smallstate = smallstate, smallchange = smallchange,
                   tstep = NA,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   newgrowthratecode = newgrowthratecode,
                   costset = costset, conjrateset, taxmattype = taxmattypeset,
                   costtype = costtype, PFrom = "PInNewSp",
                   PReplMostAbun = TRUE, duration = duration))
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

saveRDS(object = plotdata,
        file = paste0(DateTimeStamp, "_plotdata.RDS"))


#### Reading previously saved data from a CSV file ####
# # To read data from a CSV file, put the CSV file in the working directory
# # (see getwd()), uncomment this section and change the date-time-stamp in the
# # file name (the next line prints a list of files in the working directory
# # that contain 'multispecies' in their name).
# list.files(pattern = "multispecies", ignore.case = TRUE)
# # Note that the extension (.csv) should be included in the file name.
# filename <- file.path("OutputMS", "YYYY_MM_DD",
#                       "YYYY_MM_DD_MM_SSPInNewSpmultispecies.csv")
# plotdata <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
#                      dec = ".", stringsAsFactors = FALSE)
# # If plotdata has only one column, probably a semicolon instead of a comma was
# # used as separator in the CSV file. So read the file again using that separator.
# if(near(ncol(plotdata), 1L)) {
#   plotdata <- read.csv(filename, header = TRUE, sep = ";", quote = "\"",
#                        dec = ".", stringsAsFactors = FALSE)
# }
# plotdata <- as.data.frame(plotdata)
# DateTimeStamp <- substr(x = filename, start = 21, stop = 36)
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))


#### Labels and limits for plots ####
labspecies <- paste("Sp.", seq_len(maxnspecies))
names(labspecies) <- seq_len(maxnspecies)
labnspecies <- paste(nspeciesset, "sp.")
names(labnspecies) <- nspeciesset
labmodel <- c("Broken stick", "Dom. preemption")
names(labmodel) <- c(1, 2)
labcost <- paste0("Fitness cost:\n", costset, "/h")
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
plotdata$nspecies <- as.factor(plotdata$nspecies)
plotdata$intmean <- as.factor(plotdata$intmean)
plotdata$selfintmean <- as.factor(plotdata$selfintmean)
plotdata$taxmatcode <- as.factor(plotdata$taxmatcode)

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
#   theme_bw(base_size = 13) +
#   scale_x_discrete() +
#   scale_y_discrete() +
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
# 'plotdata' has not yet been converted to a dataframe, run the next line to do
# so:
# plotdata <- as.data.frame(plotdata)

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
  stop(paste("Species have different conjugation rates, so conversion",
             "from conjugationratecode\nto conjugation rate is not meaningful.",
             "\nSee the settings section in the script and in the file",
             "'settings.csv' for details"))
}
for(conjratecode_index in seq_len(conjratecode)) {
  # using dplyr::near() to allow for small (<1e-6) numeric differences
  plotdata[which(near(plotdata[, "conjratecode"], conjratecode_index)),
           "conjrate"] <- seqconjrate[conjratecode_index]
}

## Show border of ecological stability with heatmap in CreatePlot()
# Values in a facet are either all stable, or all unstable, making it impossible
# to plot contours delimiting stable and unstable regions.
CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableecol",
           filltitle = "fracstableecol", filltype = "continuous", ratio = NULL,
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode + intmean", facety = "nspecies",
           filename = "ecostabxcostyconj")


## Show border of epidemiological stability with a contour plot in CreatePlot()
# 'save' is set to FALSE and ggsave() is used to ensure the added guides
# arguments are included in the saved plots.
CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "nspecies",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode + intmean", facety = "nspecies",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmatintmeanynspecies.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "intmean",
           contour_lty = NULL,
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspecies.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}
# So intmean does not affect the border of stability in conjugation rate/cost-space

CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "nspecies",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "intmean",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmatyintmean.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "nspecies",
           contour_lty = "intmean",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = ".",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmat.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "intmean",
           contour_lty = "taxmatcode",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "selfintmean", facety = "nspecies",
           rotate_x_labels = FALSE, save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxselfintynspecies.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

##### Compare numeric thresholds with analytic predictions for two-species #####
abundance <- switch(abunmodel,
                    brokenstick = {
                      abundance <- brokenstick_fast(nspecies = 2,
                                                    totalabun = totalabun,
                                                    takelimit = TRUE)
                    },
                    dompreempt = {
                      abundance <- dompreempt_fast(nspecies = 2,
                                                   totalabun = totalabun,
                                                   takelimit = TRUE)
                    },
                    error("'abunmode' should be 'brokenstick' or 'dompreempt'"))

R1 <- abundance[1]
R2 <- abundance[2]
g22 <- conjrate_base

temp <- data.frame(conj_inter_init = rep(conjrate_base, nrow(plotdata)))
row_ind <- which(plotdata$taxmatcode == 2)
temp[row_ind, "conj_inter_init"] <- temp[row_ind, "conj_inter_init"] / 1e3

temp$p1 <- ((plotdata[, "conjrate"] * R1) - (g22 * R2))^2 +
  4 * temp[, "conj_inter_init"]^2 * R1 * R2
temp_plotdata <- plotdata
temp_plotdata$fracstableepi <- plotdata[, "conjrate"] * R1 + g22 * R2 + sqrt(temp$p1) <
  2 * plotdata[, "cost"]
rm(temp)
rm(row_ind)

temp_plotdata <- rbind(cbind(plotdata, type = "numeric"),
                       cbind(temp_plotdata, type = "analytic"))
CreatePlot(dataplot = temp_plotdata,
           xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "as.factor(type)",
           contour_lty = "taxmatcode",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Comparing epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = ".", facety = "nspecies",
           rotate_x_labels = FALSE, save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabynspecies_comparison.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}
rm(temp_plotdata)

CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "intmean",
           contour_lty = "selfintmean",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  labs(caption = NULL) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspeciescolltyinter.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
  ggsave(paste0(DateTimeStamp, "Fig20.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

# Note: assuming sets are chosen such that border of invasion is shown in the
# plot
CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "nspecies",
           contour_lty = "taxmatcode",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = ".", facety = ".",
           rotate_x_labels = FALSE, save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2) +
  annotate("text", label = c("Invasion of plasmid-\nbearing bacteria\n[WARNING: should include (self)intmean!",
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
           contour_lty = NULL, filltype = "continuous", ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
           rotate_x_labels = FALSE, filename = "epistabheatmap") +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)

CreatePlot(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableepi",
           filltitle = "fracstableepi", contour_var = NULL, contour_col = NULL,
           contour_lty = NULL, filltype = "continuous", ratio = NULL,
           title = NULL,
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode + intmean + selfintmean",
           facety = "nspecies", rotate_x_labels = TRUE, width = 9*1650/2,
           height = 2675, filename = "epistabheatmap_morefacets_v3")
