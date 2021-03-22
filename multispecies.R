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

# May RM. 2001.  Stability and complexity in model ecosystems. Princeton/Oxford:
# Princeton University Press.

# Tokeshi M. 1990. Niche apportionment or random assortment: species abundance
# patterns revisited. Journal of animal ecology 59(3):1129-1146.


#### To do ####
## General ##
# Add plasmid-bearing populations and conjugation to the model.

## geteqinfo() ##
# Currently the only the sign and complex part of the largest eigenvalue is used
# to determine the type of equilibrium. However, sometimes the largest
# eigenvalue does not have a complex part when (some of) the other eigenvalues
# do have a complex part. If this affects the equilibrium this should be taken
# into account.
# I check for repeated eigenvalues, but note that the pair a +/- bi are not
# repeated eigenvalues. See p. 133 of Edelstein-Keshet 2005 and section 10.4.3.3 from
# https://eng.libretexts.org/Bookshelves/Industrial_and_Systems_Engineering/Book%3A_Chemical_Process_Dynamics_and_Controls_(Woolf)


#### Optionally to do ####

## General ##
# I could move the 'niter' argument to the top functions, such that for 1000
# iterations, rnorm is called only once to generate 1000*nspecies growthrates,
# instead of being called 1000 times to generate nspecies growthrates.
# Or use replicate(...) from the apply-family.

# To check if enough iterations are used: run several times for niter iterations,
# if the variation in fraction of stable equilibria is too large, use more
# iterations.

# Create a matrix to store all data (all 'mydata'), then it becomes possible to
# plot species-specific growth rates and abundances.

# Create file to store default settings (e.g., sparsity = 0, intdistr =
# selfintdistr = "normal", intsd = selfintsd = 0.1) and write that to a .csv
# file with DateTimeStamp matching to other files.

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


#### Loading required libraries ####
library(deSolve)   # checkequilibrium calls ode() if showplot == TRUE
library(dplyr)     # checkequilibrium calls near()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(ggmulti)   # plotting multidimensional data

#### Settings and defining parameterspace ####

# Simulation settings
niter <- 100
saveplots <- TRUE

# Define parameter space
nspeciesset <- c(2, 4, 6)
abunmodelset <- c("brokenstick", "dompreempt")
intmeanset <- seq(from = -1.5, to = 1.5, by = 0.1)
selfintmeanset <- seq(from = -1.5, to = 1.5, by = 0.1)

# Settings for testing code
niter <- 5
saveplots <- FALSE
nspeciesset <- c(2, 4)
abunmodelset <- c("brokenstick", "dompreempt")
intmeanset <- seq(from = -1.5, to = 1.5, by = 0.25)
selfintmeanset <- seq(from = -1.5, to = 1.5, by = 0.25)

#### Functions ####

# Define the generalised Lotka-Volterra model. n is a vector of species
# abundances, growthrate is a vector of growth rates, intmat is a matrix of
# interaction coefficients where element mij gives the effect of species j on
# the growth rate of species i.
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
  S0 <- n[1:3]
  S1 <- n[4:6]
  
  dS0 <-  growthrate          * S0 * (1 - intmat %*% (S0 + S1)) - (conjmat %*% S1) * S0
  dS1 <- (growthrate - costs) * S1 * (1 - intmat %*% (S0 + S1)) + (conjmat %*% S1) * S0
  
  dn <- c(dS0, dS1)
  return(list(dn))
  })
}

# The broken stick model (= MacArthur faction model), I follow the description
# of Tokeshi (1990). Species abundances are proportional to the length of
# fragments of a stick that is broken randomly at nspecies - 1 points.
brokenstick <- function(nspecies, totalabun = 1, takelimit = FALSE) {
  stopifnot(length(nspecies) == 1, nspecies > 1,
           length(totalabun) == 1, totalabun > 0)
  if(takelimit == TRUE) {
    niter <- 1e3
    abunmat <- matrix(data = NA, nrow = niter, ncol = nspecies)
    for(iterindex in 1:niter) {
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

# The dominance preemption abundance model. The first species occupies (preempts)
# more than half of the total niche, and each subsequent species occupies more
# than half of the remainder. The last species is assigned all of the niche that
# remains, to obtain the user-defined total abundance. Over many iterations,
# each species preempts on average (0.5 + 1)/2 = 0.75 of the remainder, such
# that the model converges to the geometric series with k = 0.75 for all but the
# last species. This geometric model is used instead of the dominance preemption
# model if takelimit = TRUE.
dompreempt <- function(nspecies, totalabun = 1, takelimit = FALSE) {
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

# Create an interaction matrix for nspecies species. The fraction of sparse
# interspecies interactions can be set through 'sparsity', with sparsity = 0
# leading to a fully connected matrix, and sparsity = 1 leading to all
# off-diagonal entries equal to 0. Off-diagonal entries are interspecies
# interaction coefficients drawn from the distribution 'intdistr'. Diagonal
# entries are self-interactions drawn from the distribution 'selfintdistr'. To
# get fixed values for the interactions, choose the uniform distribution and
# provide the desired value both as the minimum and maximum of the range. The
# other arguments specify the distributions from which interaction coefficients
# are drawn. Element mij gives the effect of species j on the growth rate of
# species i.
getintmat <- function(nspecies, sparsity = 0,
                      intdistr = "normal", intmean = -1, intsd = 0.1,
                      intrange = c(-1, 0),
                      selfintdistr = "normal", selfintmean = -1, selfintsd = 0.1,
                      selfintrange = c(-1, 0)) {
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
         },
         uniform = {
           stopifnot(length(selfintrange) == 2, selfintrange[1] <= selfintrange[2])
           diag(intmat) <- runif(n = nspecies,
                                 min = selfintrange[1], max = selfintrange[2])
         },
         {
           warning("'selfintdistr' should be 'normal' or 'uniform'.")
           diag(intmat) <- NULL
         }
  )
  
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
  derivatives <- abundance*(growthrate + intmat %*% abundance)
  atequilibrium <- all(near(0, derivatives))
  if(printderivatives == TRUE) {
    print(paste("Derivatives:", paste0(round(derivatives, 4), collapse = ", ")),
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

geteqinfo <- function(abundance, intmat,
                      growthrate) {
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
  eqinfo <- c(eigvalRe, eigvalIm, eigvalReSign, eigvalImSign, eigvalRep)
  return(eqinfo)
}

# Function to create plots
CreatePlot <- function(dataplot = plotdata, xvar = "intmean", yvar = "selfintmean",
                       fillvar, filltitle, filltype = "discrete", limits = NULL, 
                       labx = "Mean interaction coefficient",
                       laby = "Mean selfinteraction coefficient",
                       mytag = NULL, addstamp = FALSE, diagional = "none",
                       facetx = "modelcode", facety = "nspecies", as.table = TRUE,
                       marginx = NULL, marginy = NULL,
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

# nspecies <- 4
# abunbrokenstick <- brokenstick(nspecies = nspecies, takelimit = TRUE)
# abundompreempt <- dompreempt(nspecies = nspecies, takelimit = TRUE)
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

# Create matrix to store data
nrowplotdata <- length(nspeciesset)*length(abunmodelset)*
  length(intmeanset)*length(selfintmeanset)
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = 11)
colnames(plotdata) <- c("niter", "nspecies", "modelcode",
                       "intmean", "selfintmean",
                       "mingrowthrate", "meangrowthrate", "maxgrowthrate",
                       "fracstable", "fracreal", "fracrep")

# Run simulations
rowindexplotdata <- 1
rowindexmydata <- 1
for(nspecies in nspeciesset) {
  for(abunmodel in abunmodelset) {
    print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
                 ": started at ", Sys.time()), quote = FALSE)
    if(abunmodel == "brokenstick") {
      abundance <- brokenstick(nspecies = nspecies, takelimit = TRUE)
      modelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt(nspecies = nspecies, takelimit = TRUE)
      modelcode <- 2
    }
    
    for(intmean in intmeanset) {
      for(selfintmean in selfintmeanset) {
        mydata <- matrix(data = NA, nrow = niter * nspecies, ncol = 15)
        for(iter in 1:niter) {
          intmat <- getintmat(nspecies = nspecies,
                              intmean = intmean, selfintmean = selfintmean)

          
          growthrate <- getgrowthrate(abundance = abundance, intmat = intmat)

          eqinfo <- geteqinfo(abundance = abundance, intmat = intmat,
                              growthrate = growthrate)
          
          mydata[(1 + nspecies*(iter - 1)):(nspecies*iter), ] <- cbind(
            niter = rep(niter, nspecies),
            nspecies = rep(nspecies, nspecies),
            abunmodel = rep(modelcode, nspecies),
            intmean = rep(intmean, nspecies),
            selfintmean = rep(selfintmean, nspecies),
            iter = rep(iter, nspecies),
            species = 1:nspecies,
            abundance = abundance,
            selfint = diag(intmat),
            growthrate = growthrate,
            eqinfo = matrix(rep(eqinfo, nspecies), nrow = nspecies, byrow = TRUE)
          )
        }

        colnames(mydata) <- c("niter", "nspecies", "abunmodel",
                              "intmean", "selfintmean",
                              "iter", "species", "abundance",
                              "selfint", "growthrate",
                              "eigvalRe", "eigvalIm",
                              "eigvalReSign", "eigvalImSign", "eigvalRep")
        fracstable <- length(which(mydata[, "eigvalRe"] < 0))/(nspecies*niter)
        fracreal <- length(which(mydata[, "eigvalImSign"] == 0))/(nspecies*niter)
        fracrep <- length(which(mydata[, "eigvalRep"] != 0))/(nspecies*niter)

        plotdata[rowindexplotdata, ] <- c(niter, nspecies, modelcode,
                                          intmean, selfintmean, min(mydata[, "growthrate"]),
                                          mean(mydata[, "growthrate"]), max(mydata[, "growthrate"]),
                                          fracstable, fracreal, fracrep)
        rowindexplotdata <- rowindexplotdata + 1
      }
    }
  }
}
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)


#### Showing and saving output ####
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
write.csv(plotdata, file = paste0(DateTimeStamp, "multispecies.csv"),
          quote = FALSE, row.names = FALSE)

## Labels and limits for plots ##
labspecies <- paste(nspeciesset, "species")
names(labspecies) <- nspeciesset
labmodel <- c("Broken stick model", "Dominance preemption model")
names(labmodel) <- c(1, 2)
mylabeller <- labeller(nspecies = labspecies, modelcode = labmodel,
                       .default = label_both)

plotdata <- as.data.frame(plotdata)
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
           facety = "nspecies", facetx = "modelcode", diagional = "minor")

CreatePlot(fillvar = "mingrowthrate", filltitle = "Minimum growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "nspecies", facetx = "modelcode", diagional = "minor")

CreatePlot(fillvar = "meangrowthrate", filltitle = "Mean growth rate",
           filltype = "binned", limits = limitsgrowthratebinned, 
           facety = "nspecies", facetx = "modelcode", diagional = "minor")

CreatePlot(fillvar = "meangrowthrate", filltitle = "Mean growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "nspecies", facetx = "modelcode", diagional = "minor")

CreatePlot(fillvar = "maxgrowthrate", filltitle = "Max growth rate",
           filltype = "binned", limits = limitsgrowthratebinned, 
           facety = "nspecies", facetx = "modelcode", diagional = "minor")

CreatePlot(fillvar = "maxgrowthrate", filltitle = "Max growth rate",
           filltype = "continuous", limits = limitsgrowthrate, 
           facety = "nspecies", facetx = "modelcode", diagional = "minor")

## Plot equilibrium characteristics
CreatePlot(fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies", facetx = "modelcode", diagional = "both")

CreatePlot(fillvar = "fracreal", filltitle = "Fraction real",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies", facetx = "modelcode", diagional = "both")

CreatePlot(fillvar = "fracrep", filltitle = "Fraction repeated eigenvalues",
           filltype = "continuous", limits = limitsfraction, 
           facety = "nspecies", facetx = "modelcode", diagional = "both")

### To test plots without using CreatePlot() ###
ggplot(data = plotdata, aes(x = intmean, y = selfintmean, fill = fracstable)) +
  geom_raster() +
  scale_x_continuous() +
  scale_y_continuous() +
  scale_fill_viridis_c(limits = limitsfraction) +
  geom_abline(intercept = 0, slope = -1, col = "white", size = 1.1) +
  geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1) +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(caption = paste(niter, "iterations")) +
  facet_grid(nspecies ~ modelcode, labeller = mylabeller)


### Show species-specific information for the last combination of intmean and
# selfintmean.
# NOTE: I use mydata because species-specific data is not in plotdata. As a
# consequence, the data shown is only for the last series of values.
mydata <- as.data.frame(mydata)

# Check relation between abundance and growthrate for the different species.
ggplot(data = mydata, aes(x = abundance, y = growthrate)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_point(aes(color = iter), size = 2) +
  scale_color_viridis_c()

# Alternatively use something like
lattice::xyplot(meangrowthrate ~ intmean | as.factor(modelcode),
                data = plotdata, groups = selfintmean)

# Using intmean instead of iter for the colorscale is more informative, but is
# not yet possible because mydata contains a single value for intmean.

ggplot(data = mydata, aes(x = selfint, y = growthrate, color = eigvalRe)) +
  geom_point() +
  scale_x_continuous() +
  scale_y_continuous() +
  scale_color_viridis_c() +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(caption = paste(niter, "iterations")) +
  facet_grid(nspecies ~ abunmodel, labeller = mylabeller)

# Becomes interesting if some eigenvalues have imaginairy part
ggplot(data = mydata, aes(x = eigvalRe, y = eigvalIm, color = selfint)) +
  geom_point() +
  scale_x_continuous() +
  scale_y_continuous() +
  scale_color_viridis_c() +
  coord_fixed(ratio = 1, expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(caption = paste(niter, "iterations")) +
  facet_grid(nspecies ~ abunmodel, labeller = mylabeller)


# Using library(ggmulti)
# Maybe clearer when using log10(abundance)
ggplot(data = mydata, aes(color = as.factor(species))) + 
  coord_serialaxes(axes.sequence = c("abundance", "selfint", "growthrate")) +
  scale_color_viridis_d() +
  geom_path()

ggplot(mydata, mapping = aes(selfint = selfint,
                             growthrate = growthrate,
                             colour = factor(species))) +
  geom_path()  + 
  coord_serialaxes()

ggplot(data = mydata, aes(x = selfint, y = growthrate, color = as.factor(species))) +
  geom_point() +
  scale_x_continuous() +
  scale_y_continuous() +
  scale_color_viridis_d() +
  theme(legend.position = "bottom") +
  labs(caption = paste(niter, "iterations")) +
  facet_grid(nspecies ~ abunmodel, labeller = mylabeller)


#### Comparing abundance models ####
for(nspecies in nspeciesset) {
  comparingabundance <- data.frame(
    nspecies = rep(nspecies, 2*nspecies),
    species = as.factor(rep(1:nspecies, 2)),
    abun = c(brokenstick(nspecies = nspecies, takelimit = TRUE),
             dompreempt(nspecies = nspecies, takelimit = TRUE)),
    model = rep(c("brokenstick", "dompreempt"), each = nspecies)
  )
  
  plot1 <- ggplot(data = comparingabundance, aes(x = species, y = abun, color = model)) +
    theme_bw() +
    scale_x_discrete(limits = factor(1:max(nspeciesset))) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "bottom") +
    labs(title = "Comparing abundance models: linear scale",
         x = "Species rank", y = "Species abundance") +
    geom_line(aes(group = model), size = 1.25) +
    geom_point(size = 2)
  print(plot1)
  if(saveplots == TRUE) {
    filename <- paste0(DateTimeStamp, "compareabun", nspecies, "species.png")
    ggsave(filename)
  }
  
  plot2 <- ggplot(data = comparingabundance, aes(x = species, y = abun, color = model)) +
    theme_bw() +
    scale_x_discrete(limits = factor(1:max(nspeciesset))) +
    scale_y_continuous(limits = c(NA, 1), trans = "log10") +
    theme(legend.position = "bottom") +
    labs(title = "Comparing abundance models: logarithmic scale",
         x = "Species rank", y = "Log10(species abundance") +
    geom_line(aes(group = model), size = 1.25) +
    geom_point(size = 2)
  print(plot2)
  if(saveplots == TRUE) {
    filename <- paste0(DateTimeStamp, "compareabun", nspecies, "specieslog.png")
    ggsave(filename)
  }
}
