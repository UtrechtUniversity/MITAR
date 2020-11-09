#### Pair-formation model of conjugation ####


#### To do ####

## Costs in growth ##
# Rethink the way of defining donor growth rate and costs: unless it is assumed 
# that the donor is the same species as the recipient, it does not make sense
# to define donor growth rates as recipient growth rate altered by costs from
# plasmid carriage. Instead, the donor could be assumed to have a growth rate
# unrelated to the recipient growth rate, and if the donor is adapted to the
# plasmid it does not have plasmid costs, but the new environment could be assumed
# to lead to costs in growth as well.

# Using (1 - c)*b to model costs implies that a plasmid decreases the growth rate
# by a certain percentage. So if you model more species and their growth rates
# are different, using the same value for c leads to different absolute decrease
# in growth rates.

# Using cd = ct is not biologically realistic.

## Other
# If I want to include stochastics, see Price 'An efficient moments-based inference
# method for within-host bacterial infection dynamics' and Volkova 'Modelling dynamics
# of plasmid-gene mediated antimicrobial resistance in enteric bacteria using stochastic
# differential equations' for ideas

# Also return ComplexEigVal and ComplexEigValBulk from function CalcEigenvalues

# Does using stol = 1.25E-6 for 8 equations in SimulationPairs and stol = 2.5E-6
# for 4 equations in SimulationBulk make sense if one wants the simulations to
# terminate at the same densities?

# For PlotOverTime, I want the plot to be always shown, and saved only if save = TRUE.
# Change to using ggplot instead of matplot.deSolve to enable saving an object through ggsave(...)
# Displaying graph first and then use if (save == TRUE) {dev.copy(png,'myplot.png'); dev.off()}
# might work without having to code the creation of the plot twice ?
# In saveplots == 1 does not work well with RunOverTime, since only the first two plots
# are saved and then warning is issued that file already exists. Add column with
# rownumber to dataframe to append to filename?

# The calculations of TotalDEstConjBulkDonor, TotalREstConjBulkDonor, TotalTransEstConjBulkDonor
# ect. for approximating bulk-rates can be moved inside the EstConjBulk function.

# Prevent overlapping values for the colorbar (could use angle with hjust, vjust)

## Voor plots over time waar ik het pair-formation en het bulk-model met elkaar Vergelijk
# is plotten van TotalD, TotalR, TotalTrans van het pair-model en DBulk, RBulk, TransBulk van het bulk-model
# een betere vergelijking, zie PlotOverTime uit het script voor het twee-compartimenten model.

# For PlotOverTime: use \n to print subtitle over two lines, use oma (see ?par and ?mtext)
# to prevent text from overlapping with legend

# If save = TRUE for plots over time, dev.off is not called?

# SummaryPlot() does not use the names of the arguments for creating titles

# aes_string is soft-deprecated (see help(aes_string)), use tidy evaluation idioms instead,
# see the quasiquotation section in aes() documentation and https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
# See also # On aes_string see https://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r

# See ggplot2::expand_limits to have the same limits and colorscale for the two plots

# Create function to filter on small negative and small positive state values and set them to 0

# In functions use more informative name arguments instead of e.g. mydf = mydf.

# Make better comparison to check for unstable equilibrium where invasion is not possible,
# currently comparison ofvalues that will not be exact

# Use dplyr to rename columns?

# Change facetting in CreatePlot() to only use facets for variables that are
# not unique. Then moving CreatePlot(fillvar = "NutrEq") and
# CreatePlot(fillvar = "log10(gtbulk)", limits = limitsbulkrates) to directly
# after they have been created will prevent plotting them multiple times because
# of different values of cd, ct, gd, gt which do not influence log10(gtbulk).
# Currently that is not possible because CreatePlot() needs cd, ct, gd, gt for
# facetting.

# Warn if lenght of variables that are not passed on to the CreatePlot function
# are > 1. E.g., if DInitset <- c(100, 1000) there will be 2 values for each
# kp*kn*cd*ct combination. I don't know how these are handled when plotting,
# probably the different values are plotted on top of each other so only the
# last value is visible. See the To-do section in the script for the
# Two-compartment model for a stub to check which parameters have multiple
# values that are not handled as such when plotting.

# Run bulk-model for short time (same as short pair-model) and compare them
# to see if the bulk-conjugation parameter holds (automated analysis somehow)

# QUESTION: have the resources to be considered for stability analysis?
# Imran does exclude them from stability analysis?
# Since I am only interested in perturbations of the equilibrium by adding donors,
# can the stability analysis be simplified (by only looking at the eigenvectors
# (eigenvalue?) of the donor-equation?). See 'A Practical Guide to Ecological
# Modelling: Using R as a Simulation Platform', p. 229-230)

# Other way of plotting data: instead of using cdSet = c(0.01, 0.025, 0.05) and
# ctSet = c(0.01, 0.025, 0.05) and plotting facets cd ~ ct, one could also use
# cdSet = 0.025 and plot facets for ct < cd, ct = cd, ct > cd. The same approach
# can be used for gd and gt.

# Rethink the comments: comment why you do something, not what you did (if you
# have to comment what you did, rewrite your code to make it clearer from the
# code itself).

# Check for use of tabs and double spaces and switch to using double spaces ?

# Look for 'hardcoded', 'hard-coded' and rethink them.

# Summary(...) is not reliable to PRINT range of data if range is > 7 orders
# of magnitude?:
# summary(c(1e-15, 1e-12, 1e-7)) prints 0.00e+00 as minimum
# summary(c(1e-15, 1e-12, 1e-7))["Min."] prints 1e-15
# summary(c(1e-15, 1e-12, 1e-8)) prints correct value

# The functions RunOverTime and PlotOverTime in the two-compartment script are
# more elaborate to enable comparison of D, R, Trans in the bulk-conjugation model
# with TotalD, TotalR, TotalTrans in the pair-formation model.
# In the two-compartment model I have not used root- and eventfunctions,
# maybe should also delete them here, since it sometimes takes very long to
# execute code (if populations go extinct and therefor many roots are found?)

## NUTRIENTS are also in the root- and event-functions, see comment at their
# function definitions

## Remarks
# Indexing of tibbles: use MyData[[1, "REq"]] to return a vector
# (MyData[1, "REq"] returns a LIST). See is.vector(state) and is.numeric(state).
# Using MyData$DInit[i] in a loop does work with tibbles as well.

## Previously names of the object returned by apply were sometimes not preserved,
# because return had a colnames attribute instead of a names attribute. Now I
# included setting names in the function, did not slow down code execution for
# 7488 rows (actually it was faster, 97 instead of 101 seconds)

#### Introduction ####
# The pair-formation model was taken from equation (2) of Zhong (2010) and
# expanded to to include costs, washout, and nutrients. See also Zhong (2012).


#### References ####

# Zhong 2010: Zhong X, Krol JE, Top EM, Krone SM. 2010. Accounting for mating
# pair formation in plasmid population dynamics. Journal of Theoretical Biology
# 262:711-719.

# Zhong 2012: Zhong X, Droesch J, Fox R, Top EM, Krone SM. 2012. On the meaning
# and estimation of plasmid transfer rates for surface-associated and well-mixed
# bacterial populations. Journal of Theoretical Biology 294:144-152.


#### Loading packages ####
library(deSolve) # Solve differential equations with results over time.
library(ggplot2) # For plotting data
library(RColorBrewer) # For better color schemes
library(rootSolve) # Integration, obtaining jacobian matrix and eigenvalues.
library(tidyr) # for 'expand.grid()' with dataframe as input
library(dplyr) # for mutate() to create new variables in dataframe/tibble

#### Plotting and simulation options ####
saveplots <- 0
atol <- 1e-10 # lower absolute error tolerance of integrator used by runsteady()
# to prevent 'DLSODE-  Warning..internal T (=R1) and H (=R2) are [1] 0 such that
# in the machine, T + H = T on the next step  [1] 0 (H = step size). Solver will
# continue anyway', which eventually leads to aborted integration.
tmaxsteady <- 1e8
timesEstConj <- seq(from = 0, to = 3, by = 0.1)
MyColorBrew <- rev(brewer.pal(11, "Spectral")) # examples: display.brewer.all()

#### Functions ####
# Calculate the plasmid-free equilibrium (R*, Nutr*) using the solution to
# dR/dt = R*(Nutr*bR - w) == 0, dNutr/dt = (NI - Nutr)*w - NutrConv*Nutr*bR*R ==
# 0 with R > 0 and Nutr > 0
CalcEqPlasmidfree <- function(MyData) {
  with(as.list(MyData), {
    REq <- (NI - (w / bR)) / NutrConv
    NutrEq <- w / bR
    Eq <- c(NutrEq = NutrEq, REq = REq)
    return(Eq)
  })
}

# ODE-model used to approximate the bulk-conjugation rate of the donor.
# Nutrients, growth, washout, conjugation from transconjugants, and Mtt-pairs
# are not included in this model.
ModelEstConjBulkDonor <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dD <- - kp*D*R + kn*(Mdr + Mdt)
    dR <- - kp*R*(D + Trans) + kn*(Mdr + Mrt)
    dTrans <- - kp*R*Trans + kn*(Mdt + Mrt)
    dMdr <- kp*D*R - kn*Mdr - gd*Mdr
    dMdt <- gd*Mdr - kn*Mdt
    dMrt <- kp*R*Trans - kn*Mrt
    return(list(c(dD, dR, dTrans, dMdr, dMdt, dMrt)))
  })
}

# ODE-model used to approximate the bulk-conjugation rate of the transconjugant.
# Nutrients, growth, washout, and donors are not included in this model.
ModelEstConjBulkTrans <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR <- - kp*R*Trans + kn*Mrt
    dTrans <- - kp*R*Trans + kn*(Mrt + 2*Mtt)
    dMrt <- kp*R*Trans - kn*Mrt - gt*Mrt
    dMtt <- gt*Mrt - kn*Mtt
    return(list(c(dR, dTrans, dMrt, dMtt)))
  })
}

# Function to estimate bulk-conjugation rates by running simulations with the
# adjusted pair-formation models for a short time (i.e., tail(timesEstConj, 1)
# hours) and calculate approximations of gdbulk and gtbulk from the output,
# following Zhong's approach for the calculations.
EstConjBulk <- function(MyData) {
  with(as.list(MyData), {
    state <- c(D = MyData[["DInit"]], R = MyData[["REq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)

    state <- c(R = MyData[["REq"]], Trans = MyData[["DInit"]], Mrt = 0, Mtt = 0)
    DataEstConjBulkTrans <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)

    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
  })
}

# ODE-model describing pair-formation and conjugation. Pair-formation between
# plasmid-free recipients and plasmid-bearing donors or transconjugants depends
# on attachment rate kp. Conjugation from donors or transconjugants occurs in
# the Mdr and Mrt pairs with intrinsic conjugation rates gd and gt respectively.
# This leads to formation of Mdt and Mtt pairs. Pairs fall apart with detachment
# rate kn. This structure of pair-formation is based on Zhong's model (Zhong
# 2010). I expanded the model to include costs in growth for plasmid-bearing
# bacteria, washout, and nutrients.
ModelPairsNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w - NutrConv*Nutr*((1 - cd)*bR*(D + Mdr + Mdt) + bR*(R + Mdr + Mrt) +
                                              (1 - ct)*bR*(Trans + Mdt + Mrt + 2*Mtt))
    dD <- (1 - cd)*bR*Nutr*(D + Mdr + Mdt) - kp*D*R + kn*(Mdr + Mdt) - w*D
    dR <- bR*Nutr*(R + Mdr + Mrt) - kp*R*(D + Trans) + kn*(Mdr + Mrt) - w*R
    dTrans <- (1 - ct)*bR*Nutr*(Trans + Mdt + Mrt + 2*Mtt) - kp*R*Trans + kn*(Mdt + Mrt + 2*Mtt) -
      w*Trans
    dMdr <- kp*D*R - kn*Mdr - gd*Mdr - w*Mdr
    dMdt <- gd*Mdr - kn*Mdt - w*Mdt
    dMrt <- kp*R*Trans - kn*Mrt - gt*Mrt - w*Mrt
    dMtt <- gt*Mrt - kn*Mtt - w*Mtt
    return(list(c(dNutr, dD, dR, dTrans, dMdr, dMdt, dMrt, dMtt)))
  })
}

# The bulk-conjugation model
ModelBulkNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w - NutrConv*bR*Nutr*((1 - cd)*D + R + (1 - ct)*Trans)
    dD <- ((1 - cd)*bR*Nutr - w)*D
    dR <- (bR*Nutr - gdbulk*D - gtbulk*Trans - w)*R
    dTrans <- (1 - ct)*bR*Nutr*Trans + gdbulk*D*R + gtbulk*Trans*R - w*Trans
    return(list(c(dNutr, dD, dR, dTrans)))
  })
}

# Numerically estimate the Jacobian matrix of the plasmid-free equilibrium of
# the models, then calculate (or approximate?) the eigenvalues of this matrix.
# The maximum real part of the eigenvalues is used to determine stability.
CalcEigenvalues <- function(MyData) {
  parms <- MyData
  EqFull <- c(Nutr = MyData[["NutrEq"]], D = 0, R = MyData[["REq"]], Trans = 0,
              Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  EigVal <- eigen(x = jacobian.full(y = EqFull, func = ModelPairsNutr, parms = parms),
                  symmetric = FALSE, only.values = TRUE)$values
  ComplexEigVal <- is.complex(EigVal) 
  EigVal <- Re(EigVal)
  names(EigVal) <- paste0("Eigval", 1:length(EigVal))
  DomEigVal <- max(EigVal)
  SignDomEigVal <- sign(DomEigVal)
  SignEigValEqual <- identical(rep(SignDomEigVal, length(EigVal)), sign(Re(EigVal)))
  
  EqFullBulk <- c(Nutr = MyData[["NutrEq"]], D = 0, R = MyData[["REq"]], Trans = 0)
  EigValBulk <- eigen(x = jacobian.full(y = EqFullBulk, func = ModelBulkNutr, parms = parms),
                      symmetric = FALSE, only.values = TRUE)$values
  ComplexEigValBulk <- is.complex(EigValBulk) 
  EigValBulk <- Re(EigValBulk)
  names(EigValBulk) <- paste0("EigvalBulk", 1:length(EigValBulk))
  DomEigValBulk <- max(EigValBulk)
  SignDomEigValBulk <- sign(DomEigValBulk)
  SignEigValEqualBulk <- identical(rep(sign(DomEigValBulk), length(EigValBulk)),
                                   sign(Re(EigValBulk)))
  InfoEigVal <- c(EigVal, DomEigVal = DomEigVal, SignDomEigVal = SignDomEigVal,
                  SignEigValEqual = SignEigValEqual,
                  EigValBulk, DomEigValBulk = DomEigValBulk,
                  SignDomEigValBulk = SignDomEigValBulk,
                  SignEigValEqualBulk = SignEigValEqualBulk)
  return(InfoEigVal)
}

# Simulations using the pair-formation and the bulk model. The plasmid-free
# equilibrium (Nutr*, R*) with the addition of DInit donor bacteria per mL is
# used as initial state. Note that stol is based on the average of absolute
# rates of change, not the sum.
SimulationPairs <- function(InputSimulationPairs) {
  parms <- InputSimulationPairs
  state <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]],
             R = parms[["REq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelPairsNutr,
                   parms = parms, stol = 1.25e-6, atol = atol)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

SimulationBulk <- function(InputSimulationBulk) {
  parms <- InputSimulationBulk
  state <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]],
             R = parms[["REq"]], Trans = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelBulkNutr,
                   parms = parms, stol = 2.5e-6, atol = atol)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

# Create heatmaps, save if needed
CreatePlot <- function(fillvar, gradient2 = 0, limits = NULL, midpoint = 0, dataplot = MyData,
                       xvar = "log10(kp)", yvar = "log10(kn)",
                       facetx = "cd + gd", facety = "ct + gt", save = saveplots, ...) {
  if(exists("DateTimeStamp") == FALSE) {
    warning("DateTimeStamp created to include in plot but does not correspond to filename of the dataset")
    DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
  }
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar)) + 
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(as.formula(paste(facetx, "~", facety)), labeller = label_both) +
    labs(caption = DateTimeStamp) +
    theme(legend.position = "bottom", plot.caption = element_text(vjust = 20))
  if(gradient2 == 1) {
    p <- p + scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = midpoint, limits = limits)
  } else {
    p <- p + scale_fill_gradientn(colours = MyColorBrew, limits = limits)
    # p <- p + scale_fill_distiller(palette = "Spectral", direction = 1, limits = limits)
  }
  print(p)
  if(save == TRUE) {
    fillvarname <- gsub("/", ".", fillvar)
    fillvarname <- gsub(" ", "", fillvarname)
    filename <- paste0(DateTimeStamp, "output", fillvarname, ".png")
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      ggsave(filename)
    }
  }
}

# Summarize the output per variable, indicate the orders of magnitude difference
# and negative values
SummaryPlot <- function(plotvar = plotvar, sortvalues = FALSE, ylim = NULL) {
  print(summary(plotvar))
  print("Range:", quote = FALSE)
  print(range(plotvar))
  plotdf <- data.frame(val = plotvar, sign = sign(plotvar), absval = abs(plotvar), plotcol = "black")
  if(sortvalues == TRUE) {
    plotdf <- plotdf[order(plotdf$sign, plotdf$absval), ]
  }
  valueszero <- which(plotdf[, "val"] == 0)
  nvalueszero <- length(valueszero)
  
  if(nvalueszero > 0) {
    plotdf[valueszero, "plotcol"] <- "blue"
    
    if(nvalueszero == length(plotvar)) {
      print("All values are 0!", quote = FALSE)
      minval <- 1
      maxval <- 1
    } else {
      minabsval <- min(plotdf[-valueszero, "absval"])
      maxabsval <- max(plotdf[-valueszero, "absval"])
      print(paste("After removing the", nvalueszero, "values equal to zero:"), quote = FALSE)
      print(summary(plotdf[-valueszero, "val"]))
    }
  } else {
    minabsval <- min(plotdf[, "absval"])
    maxabsval <- max(plotdf[, "absval"])
  }
  
  valuesnegative <- which(plotdf[, "val"] < 0)
  nvaluesnegative <- length(valuesnegative)
  if(nvaluesnegative > 0) {
    plotdf[valuesnegative, "plotcol"] <- "red"
    print(paste("Summary of the", nvaluesnegative, "smaller than zero:"), quote = FALSE)
    print(summary(plotdf[valuesnegative, "val"]))
    print(paste("Range of the", nvaluesnegative, "smaller than zero:"), quote = FALSE)
    print(range(plotdf[valuesnegative, "val"]))
  }
  
  ordersdiff <- log10(maxabsval / minabsval)
  
  if(ordersdiff < 2) {
    plot(plotdf[, "val"], ylim = ylim, pch = 19, col = plotdf[, "plotcol"])
  } else {
    print(paste("Using logscale because values differ over", round(ordersdiff), "orders of magnitude"), quote = FALSE)
    if(nvalueszero > 0) {
      print(paste("Not showing", nvalueszero, "values that are equal to 0"), quote = FALSE)
      plot(plotdf[-valueszero, "absval"], ylim = ylim, log = "y", pch = 19, col = plotdf[-valueszero, "plotcol"])
    } else {
      plot(plotdf[, "absval"], ylim = ylim, log = "y", pch = 19, col = plotdf[, "plotcol"])
    }
    if(nvaluesnegative > 0) {
      print(paste("Taking the absolute value of the", nvaluesnegative, "values smaller than zero"), quote = FALSE)
    }
  }
  grid()
  abline(h = 1)
  abline(h = 0.001)
}


#### Parameter values ####

## Explanation of parameters
# Growthrates (1/h) for recipient: bR (growth rates for donor and transconjugant
# are not explicitly defined, but calculated as (1 - cd)*bR and (1 - ct)*bR)
# Conjugation rates (1/h) from donor and transconjugant: gd and gt
# Mating pair attachment rate (mL * cell^-1 * h^-1): kp (k+ in Zhong's notation)
# Mating pair detachment rate (1/h): kn (k- in Zhong's notation)
# Nutrient concentration in the inflowing liquid (microgram * mL^-1): NI
# Resource conversion rate (microgram per cell division): e
# Washout rate (1/h): w

## To read data from csv-file
# FileName <- "2020_augustus_03_12_54_32outputnosimulationtwocompartment.csv"
# MyData <- read.csv(FileName, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE)
# MyData <- as.data.frame(MyData)
# DateTimeStamp <- substr(FileName, 1, nchar(FileName) - 36)

## Stable co-existence of recipients, transconjugants, and Mrt and Mtt pairs
# bRSet <- 1.7; NISet <- 10; kpSet <- 10^-9.6; knSet <- 10^0.5
# NutrConvSet <- 1E-6; wSet <- 0.04; cdSet <- 0.05; ctSet <- 0.05;
# gdSet <- 10^1.176; gtSet <- 10^1.176; DInitSet <- 1000
# Leads to time = 843129.4, Nutr = 0.02427729 D = 0, R = 3829217, Trans = 6142815,
# Mdr = Mdt = 0, Mrt = 324.6586, Mtt = 1520.435, timeBulk = 728451.6,
# NutrBulk = 0.02430818, DBulk = 0, RBulk = 3583808, TransBulk = 6391884.

## Parameterset 1: show NI and w influence stability of plasmid-free equilibrium
DInitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- 10^seq(from = 0, to = 2, by = 1)
NutrConvSet <- 1e-6
# ln(2)/24, 1/24, 1% remaining after 24h, 0.1% remaining after 24h:
wSet <- c(0.029, 0.042, 0.192, 0.288) 
kpSet <- 10^seq(from = -12, to = -6, by = 0.25)
knSet <- 10^seq(from = -1, to = 3, by = 0.25)
cdSet <- c(0.05)
ctSet <- c(0.01)
gdSet <- c(15)
gtSet <- c(15)

## Parameterset 2: show influence of costs, conjugation, attachment, and 
## detachment rates on the stability of the plasmid-free equilibrium and on
## bulk-conjugation rates.
DInitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
NutrConvSet <- 1e-6
wSet <- c(0.04)
kpSet <- 10^seq(from = -12, to = -6, by = 0.25)
knSet <- 10^seq(from = -1, to = 3, by = 0.25)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(1, 15)
gtSet <- c(1, 15)

## Using parameterset 2 with steps of 0.2 for kp and kn does not work for the
# pair-formation model, or for the bulk-conjugation model. Lowering atol to 1e-11
# doesn't resolve this. Maybe retry using jactype = "sparse", or use stode(s?)
# and/or supply jacobian if integration leads to errors ?

## Parameterset 2b: what happens at very low detachment rates?
DInitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- 10
NutrConvSet <- 1e-6
wSet <- 0.04
kpSet <- 10^seq(from = -13, to = -4, by = 0.2)
knSet <- 10^seq(from = -4, to = 3, by = 0.2)
cdSet <- c(0.05)
ctSet <- c(0.01)
gdSet <- c(15)
gtSet <- c(15)

#### Main script ####

CheckParms <- c(DInitSet, bRSet, NISet, NutrConvSet, wSet, kpSet, knSet, cdSet, ctSet)
if(any(CheckParms <= 0)) warning("All parameters should have positive values.")
if(any(c(cdSet, ctSet) >= 1)) warning("Costs should be larger than 0 and smaller than 1.")

TotalIterations <- length(DInitSet)*length(bRSet)*length(NISet)*length(NutrConvSet)*
  length(wSet)*length(kpSet)*length(knSet)*length(cdSet)*length(ctSet)*
  length(gdSet)*length(gtSet)
TotalIterations

## Calculate plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(bR = bRSet, NI = NISet, NutrConv = NutrConvSet, w = wSet)

Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree))
MyData <- cbind(MyData, Eqplasmidfree)

if(any(Eqplasmidfree[, "REq"] <= 0)) {
  warning("The number of recipients at equilibrium is not always positive.
  Rows with negative densities have been discarded!
  Increase growth rate or the nutrient concentration in the inflowing liquid, or decrease the outflow rate to prevent this!")
  MyData <- MyData[-which(MyData[, "REq"] < 0), ]
}

print("Plasmid-free equilibrium calculated:")
print(Sys.time())

## Add combinations with the parameters needed to approximate gdbulk and gtbulk
MyData <- expand_grid(MyData, gd = gdSet, gt = gtSet, DInit = DInitSet,
                      kp = kpSet, kn = knSet)
DataEstConjBulk <- t(apply(X = MyData, MARGIN = 1, FUN = EstConjBulk))

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulk <- unname(MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor))

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] + DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] + DataEstConjBulk[, "TransMrt"]
gtbulk <- unname(MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans))

MyData <- cbind(MyData, gdbulk = gdbulk, gtbulk = gtbulk)
print("Bulk-conjugation rates estimated:")
print(Sys.time())

# Approximate eigenvalues
MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)

MyInfoEigVal <- t(apply(MyData, MARGIN = 1, FUN = CalcEigenvalues))
MyData <- cbind(MyData, MyInfoEigVal)
print("Eigenvalues estimated:")
print(Sys.time())

DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
write.csv(MyData, file = paste0(DateTimeStamp, "outputnosimulation.csv"),
          quote = FALSE, row.names = FALSE)


#### Plotting output for parameterset 1 ####

# Show influence of washout rate and inflowing nutrient concentration on
# stability of the plasmid-free equilibrium (run with parameterset 1).
# Note: hardcoded legend and axis labels
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal))) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(w ~ NI, labeller = label_both) +
  labs(caption = DateTimeStamp, x = "log10(attachment rate)",
       y = "log10(detachment rate)") +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("1" = "darkred", "-1" = "darkblue"),
                    name = "Plasmid can invade",
                    labels = c("No", "Yes"))
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "outputfactor(SignDomEigValNIw).png"))
}

# These plots (not shown in article) show that nutrient inflow influences
# recipient cell density, whereas washout rate influences nutrient concentration
# (run with parameterset 1).
CreatePlot(fillvar = "log10(REq)", facetx = "w", facety = "NI")
CreatePlot(fillvar = "NutrEq", facetx = "w", facety = "NI")

# Show that dominant eigenvalues have equal signs for pair-formation and bulk
# model (run with parameterset 1, plot not shown in article)
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal == SignDomEigValBulk))) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(w ~ NI, labeller = label_both) +
  labs(caption = DateTimeStamp, x = "log10(attachment rate)",
       y = "log10(detachment rate)") +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("TRUE" = "aquamarine", "FALSE" = "red"),
                    name = "Dominant eigenvalues \nhave equal signs")
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "DifferenceInSignEigenvaluesNIw.png"))
}


#### Plotting output for parameterset 2 ####

# To show influence of costs and intrinsic conjugation rates on stability of the
# plasmid-free equilibrium, run with multiple values for kn, kp, cd, ct, gd, gt
# and than plot:
# Note: hardcoded legend
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal))) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
  labs(caption = DateTimeStamp, x = "log10(attachment rate)",
       y = "log10(detachment rate)") +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("1" = "darkred", "-1" = "darkblue"),
                    name = "Plasmid can invade",
                    labels = c("No", "Yes"))
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "outputfactor(SignDomEigVal).png"))
}

# Stability of the equilibrium for the bulk-conjugation model
# (plot not shown in article)
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigValBulk))) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
  labs(caption = DateTimeStamp, x = "log10(attachment rate)",
       y = "log10(detachment rate)") +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("1" = "darkred", "-1" = "darkblue"),
                    name = "Plasmid can invade",
                    labels = c("No", "Yes"))
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "outputfactor(SignDomEigValBulk).png"))
}

# Show if sign of dominant eigenvalues for pair-formation and bulk model is the same 
# (Figure S4 in article)
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal == SignDomEigValBulk))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
  labs(caption = DateTimeStamp) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("TRUE" = "aquamarine", "FALSE" = "red"),
                    name = "Dominant eigenvalues \nhave equal signs")
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "DifferenceInSignEigenvalues.png"))
}

## The influence of conjugation, attachment, and detachment rates on the 
## bulk-conjugation rates. Data is filtered to show only one value for costs,
## because costs do not influence the bulk-conjugation rates (Figure 4 in the
## article).
limitsbulkrates <- c(floor(min(log10(c(MyData$gdbulk, MyData$gtbulk)))),
                     ceiling(max(log10(c(MyData$gdbulk, MyData$gtbulk)))))
CreatePlot(dataplot = filter(MyData, gt == 15 & cd == 0.05 & ct == 0.01),
           fillvar = "log10(gdbulk)", facetx = "gt", facety = "gd",
           limits = limitsbulkrates)
CreatePlot(dataplot = filter(MyData, gd == 15 & cd == 0.05 & ct == 0.01),
           fillvar = "log10(gtbulk)", facetx = "gd", facety = "gt",
           limits = limitsbulkrates)

CreatePlot(fillvar = "log10(gdbulk)", facetx = ".", facety = ".",
           limits = limitsbulkrates, save = FALSE)

print("Finished plotting:")
print(Sys.time())

# If invasion is possible, run simulation to see how many bacteria of each
# population are present at equilibrium
IndexSimulation <- which(MyData$SignDomEigVal != -1)
print(paste(length(IndexSimulation), "simulations to run for the pair-formation model"))
ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))
InputSimulationPairs <- MyData[IndexSimulation, ColumnsToSelect]
OutputSimulationPairs <- t(apply(X = InputSimulationPairs, MARGIN = 1,
                                 FUN = SimulationPairs))

if(length(IndexSimulation) < nrow(MyData)) {
  NoSimulationNeeded <- cbind(time = 0, steady = 1,
                              Nutr = MyData[-IndexSimulation, "NutrEq"],
                              D = 0, R = MyData[-IndexSimulation, "REq"],
                              Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  MyData <- rbind(cbind(MyData[IndexSimulation, ], OutputSimulationPairs),
                  cbind(MyData[-IndexSimulation, ], NoSimulationNeeded))
} else {
  MyData <- cbind(MyData, OutputSimulationPairs)
}
if(any(MyData$steady == 0)) warning("Steady-state has not always been reached")

print("Pair-formation model completed running:")
print(Sys.time())

IndexSimulationBulk <- which(MyData$SignDomEigValBulk != -1)
print(paste(length(IndexSimulationBulk), "simulations to run for the bulk model"))
InputSimulationBulk <- MyData[IndexSimulationBulk, ColumnsToSelect]
OutputSimulationBulk <- t(apply(X = InputSimulationBulk, MARGIN = 1, FUN = SimulationBulk))
colnames(OutputSimulationBulk) <- paste0(colnames(OutputSimulationBulk), "Bulk")

if(length(IndexSimulationBulk) < nrow(MyData)) {
  NoSimulationNeededBulk <- cbind(timeBulk = 0, steadyBulk = 1, NutrBulk = MyData[-IndexSimulationBulk, "NutrEq"],
                                  DBulk = 0, RBulk = MyData[-IndexSimulationBulk, "REq"],
                                  TransBulk = 0)
  MyData <- rbind(cbind(MyData[IndexSimulationBulk, ], OutputSimulationBulk),
                  cbind(MyData[-IndexSimulationBulk, ], NoSimulationNeededBulk))
} else {
  MyData <- cbind(MyData, OutputSimulationBulk)
}
if(any(MyData$steadyBulk == 0)) warning("Steady-state has not always been reached")

print("Bulk-conjugation model completed running:")
print(Sys.time())

MyData <- cbind(MyData, TotalD = NA, TotalR = NA, TotalTrans = NA, TotalPlasmid = NA, TotalBio = NA,
                TotalPlasmidBulk = NA, TotalBioBulk = NA)
MyData[, "TotalD"] <- MyData[, "D"] + MyData[, "Mdr"] + MyData[, "Mdt"]
MyData[, "TotalR"] <- MyData[, "R"] + MyData[, "Mdr"] + MyData[, "Mrt"]
MyData[, "TotalTrans"] <- MyData[, "Trans"] + MyData[, "Mdt"] + MyData[, "Mrt"] + 2*MyData[, "Mtt"]
MyData[, "TotalPlasmid"] <- MyData[, "TotalD"] + MyData[, "TotalTrans"]
MyData[, "TotalBio"] <- MyData[, "TotalR"] + MyData[, "TotalPlasmid"]
MyData[, "TotalPlasmidBulk"] <- MyData[, "DBulk"] + MyData[, "TransBulk"]
MyData[, "TotalBioBulk"] <- MyData[, "RBulk"] + MyData[, "TotalPlasmidBulk"]

write.csv(MyData, file = paste0(DateTimeStamp, "outputsimulation.csv"),
          quote = FALSE, row.names = FALSE)

#### Plotting summaries of the variables and of some parameters ####
SummaryPlot(MyData$NutrEq)
SummaryPlot(MyData$REq)
SummaryPlot(MyData$gdbulk)
SummaryPlot(MyData$gtbulk)
SummaryPlot(MyData$DomEigVal)
SummaryPlot(MyData$DomEigValBulk)
SummaryPlot(MyData$Nutr)
SummaryPlot(MyData$D)
SummaryPlot(MyData$R)
SummaryPlot(MyData$Trans)
SummaryPlot(MyData$Mdr)
SummaryPlot(MyData$Mdt)
SummaryPlot(MyData$Mrt)
SummaryPlot(MyData$Mtt)
SummaryPlot(MyData$NutrBulk)
SummaryPlot(MyData$DBulk)
SummaryPlot(MyData$RBulk)
SummaryPlot(MyData$TransBulk)
SummaryPlot(MyData$TotalD)
SummaryPlot(MyData$TotalR)
SummaryPlot(MyData$TotalTrans)
SummaryPlot(MyData$TotalPlasmid) # Plots on linear scale despite varying over many orders of magnitude ?
SummaryPlot(MyData$TotalPlasmidBulk)

### Some controls
# Equilibrium at t=0 despite unstable plasmid-free equilibrium
if(any(MyData$time == 0 & MyData$SignDomEigVal == 1)) {
  A <- length(which(MyData$time == 0 & MyData$SignDomEigVal == 1))
  warning(paste("Equilibrium at t=0 despite unstable plasmid-free equilibrium in", A, "cases!"))
}

if(any(MyData$timeBulk == 0 & MyData$SignDomEigVal == 1)) {
  A <- length(which(MyData$timeBulk == 0 & MyData$SignDomEigVal == 1))
  warning(paste("Equilibrium at t=0 despite unstable plasmid-free equilibrium in", A, "cases!"))
}

# Unstable equilibrium, but invasion not possible
if(any(MyData$TotalR == MyData$TotalBio & MyData$SignDomEigVal == 1)) {
  A <- length(which(MyData$TotalR == MyData$TotalBio & MyData$SignDomEigVal == 1))
  warning(paste("Unstable equilibrium but invasion not possible with pair-model in", A, "cases!"))
}

if(any(MyData$RBulk == MyData$TotalBioBulk & MyData$SignDomEigValBulk == 1)) {
  A <- length(which(MyData$RBulk == MyData$TotalBioBulk & MyData$SignDomEigValBulk == 1))
  warning(paste("Unstable equilibrium but invasion not possible with bulk-model in", A, "cases!"))
}

# If this differs, comparisons between the 2 models should use fractions, not cell counts
summary(MyData$TotalBio / MyData$TotalBioBulk)
range(MyData$TotalBio / MyData$TotalBioBulk)

# If the biomass is the same across simulations, one could use the number of cells instead of fractions of total biomass
summary(MyData$TotalBio)
max(MyData$TotalBio) / min(MyData$TotalBio)
summary(MyData$TotalBioBulk)
max(MyData$TotalBioBulk) / min(MyData$TotalBioBulk)

# Which proportion of cells is plasmid-bearing, if it is not nearly 0 or nearly 1
PossibleCoexistence <- MyData[which(MyData[, "TotalPlasmid"]/MyData[, "TotalBio"] > 1e-3 &
                                          MyData[, "TotalPlasmid"]/MyData[, "TotalBio"] < 0.999), ]
dim(PossibleCoexistence)
write.csv(PossibleCoexistence, file = paste0(DateTimeStamp, "possiblecoexistence.csv"),
          quote = FALSE, row.names = FALSE)

PossibleCoexistenceBulk <- MyData[which(MyData[, "TotalPlasmidBulk"]/MyData[, "TotalBioBulk"] > 1e-3 &
                    MyData[, "TotalPlasmidBulk"]/MyData[, "TotalBioBulk"] < 0.999), ]
dim(PossibleCoexistenceBulk)
write.csv(PossibleCoexistenceBulk, file = paste0(DateTimeStamp, "possiblecoexistencebulk.csv"),
          quote = FALSE, row.names = FALSE)

#### Plots for parameterset 1 ####
CreatePlot(fillvar = "TotalPlasmid / TotalBio", facetx = "w", facety = "NI")
CreatePlot(fillvar = "TotalPlasmidBulk / TotalBioBulk", facetx = "w", facety = "NI")

#### Plots for parameterset 2 ####

# Compare biomass across the two models
CreatePlot(fillvar = "TotalBio / TotalBioBulk")

# Fraction plasmid-bearing cells
CreatePlot(fillvar = "TotalPlasmid/TotalBio") # Figure S3 with parameterset 2
CreatePlot(fillvar = "TotalPlasmidBulk/TotalBioBulk")

# Fraction donors
CreatePlot(fillvar = "TotalD/TotalBio") # Figure S2 with parameterset 2
CreatePlot(fillvar = "DBulk/TotalBioBulk")

# Fraction recipients
CreatePlot(fillvar = "TotalR/TotalBio")
CreatePlot(fillvar = "RBulk/TotalBioBulk")

# Fraction transconjugants
CreatePlot(fillvar = "TotalTrans/TotalBio") # Figure S1 with parameterset 2
CreatePlot(fillvar = "TransBulk/TotalBioBulk")

################################################################################
BackupMyData <- MyData

###### To create plots over time

# Settings for simulations, plotting, and printing
mylty <- c(lty = c(3, 1, 2, 1, 1, 1, 1, 1))
# mycol <- c("black", "purple", "green1", "red", "yellow", "hotpink", "blue", "cyan")
mycol <- c("black", brewer.pal(7, "Set1"))
myylim <- c(1E-7, 1E7) # Defining the limits for the y-axis
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
plotoutput <- 1
verbose <- 0 # if verbose == 1, diagnostics on the simulations are printed and roots are indicated in the graphs
smallchange <- c(1E-5)
Mytmax <- c(1E5)
Mytstep <- c(10)

#### Create matrix to store data ####
TheseRows <- c(1, nrow(MyData))
TheseRows <- 1:nrow(MyData)
TheseRows <- ceiling(seq(from = nrow(MyData)/(4*4*2*2), to = nrow(MyData), length.out = 4*4*2))
if(!exists("ColumnsToSelect")) {ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))}
Mydf <- MyData[TheseRows, ColumnsToSelect]
TotalIterations <- length(TheseRows)
print(TotalIterations)
CurrentIteration <- 0

# Times for which output of the simulation is wanted.
# Note that the used ode-solvers are variable-step methods, so the times in times
# are NOT the only times at which integration is performed. See
# help(diagnostics.deSolve()) and help(lsodar()) for details.
times <- c(0:100, seq(from = 100 + Mytstep, to = Mytmax, by = Mytstep))

# Define root-function
# The root-argument becomes equal to 0 if the sum of absolute rates of change is
# equal to the threshold smallchange, i.e., when equilibrium is nearly reached.
# Then the simulation is terminated. See help(events) and help(lsodar) (both in
# the deSolve package) for background information and examples.
rootfun <- function(t, state, parms) {
  sum(abs(unlist(ModelPairsNutr(t, state, parms)))) - smallchange
}

rootfunBulk <- function(t, stateBulk, parmsBulk) {
  sum(abs(unlist(ModelBulkNutr(t, stateBulk, parmsBulk)))) - smallchange
}

# Note: the functions RunOverTime and PlotOverTime in the two-compartment script are
# more elaborate to enable comparison of D, R, Trans in the bulk-conjugation model
# with TotalD, TotalR, TotalTrans in the pair-formation model.
RunOverTime <- function(parms = Mydf, verbose = FALSE, ...) {
  state <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]], R = parms[["REq"]],
             Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  out2 <- ode(t = times, y = state, func = ModelPairsNutr, parms = parms,
           rootfun = rootfun, verbose = verbose)
  EqAfterInvasion <- tail(out2, 1)
  print(EqAfterInvasion)
  if(verbose == TRUE) {
    print(diagnostics(out2))
    print(attributes(out2))
  }
  PlotOverTime(plotdata = out2, parms = parms, type = "Pair", verbose = verbose, saveplot = saveplots)
  stateBulk <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]], R = parms[["REq"]], Trans = 0)
  out2bulk <- ode(t = times, y = stateBulk, func = ModelBulkNutr, parms = parms,
                rootfun = rootfunBulk, verbose = verbose)
  EqAfterInvasionBulk <- tail(out2bulk, 1)
  if(verbose == TRUE) {
    print(diagnostics(out2bulk))
    print(attributes(out2bulk))
  }
  PlotOverTime(plotdata = out2bulk, parms = parms, type = "Bulk", verbose = verbose, saveplot = saveplots)
  
  EqAfterInvasionTotal <- cbind(EqAfterInvasion, EqAfterInvasionBulk)
  names(EqAfterInvasionTotal) <- c(colnames(EqAfterInvasion),
                                   paste0(colnames(EqAfterInvasionBulk), "Bulk"))
  return(EqAfterInvasionTotal)
}

PlotOverTime <- function(plotdata = out2, parms = parms, type = "Pair", verbose = FALSE, saveplot = saveplots, ...) {
  maintitle <- paste(type, "model")
  subtitle <- paste0("bR=", parms[["bR"]], " NI=", parms[["NI"]],
                     " log10kp=", log10(parms[["kp"]]), " log10kn=", log10(parms[["kn"]]),
                     " NutrConv=", parms[["NutrConv"]], " w=", parms[["w"]],
                     " cd=", parms[["cd"]], " ct=", parms[["ct"]])
  if(type == "Pair") {
    subtitle <- paste0(subtitle, " gd=", parms[["gd"]], " gt=", parms[["gt"]]) 
  } else {
    subtitle <- paste0(subtitle, " gdbulk=",  signif(parms[["gdbulk"]], digits = 4),
                       " gtbulk=", signif(parms[["gtbulk"]], digits = 4))
  }
  if(saveplot == TRUE) {
    filename <- paste0(DateTimeStamp, "output", gsub(" ", "", maintitle), ".png")
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      png(filename = filename)
    }
    matplot.deSolve(plotdata, main = maintitle, sub = subtitle, ylim = myylim,
                    log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                    legend = list(x = "topright"))
    grid()
    if(verbose == TRUE) {
      abline(v = attributes(plotdata)$troot)
    }
    if(file.exists(filename) == FALSE) {
      dev.off()
    }
  } else {
    matplot.deSolve(plotdata, main = maintitle, sub = subtitle, ylim = myylim,
                    log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                    legend = list(x = "bottomright"))
    grid()
    if(verbose == TRUE) {
      abline(v = attributes(plotdata)$troot)
    }
  }
}

EqAfterInvasion <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime))

EqAfterInvasion <- cbind(Mydf, EqAfterInvasion)

write.csv(EqAfterInvasion, file = paste0(DateTimeStamp, "outputrunovertimeonecomp.csv"),
          quote = FALSE, row.names = FALSE)

#### The code below is not used ####
# MyData <- matrix(data = NA, nrow = TotalIterations, ncol = 61, byrow = TRUE) # To store output data

EqAfterInvasionTotal <- cbind(EqAfterInvasionTotal,
                              TotalTrans = EqAfterInvasionTotal[, "Trans"] +
                              EqAfterInvasionTotal[, "Mdt"] +
                              EqAfterInvasionTotal[, "Mrt"] +
                              2*EqAfterInvasionTotal[, "Mtt"])

EqAfterInvasionTotal <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime))


MyData[CurrentIteration, c(1:14)] <- unname(c(parmsBulk, Eq))
MyData[CurrentIteration, c(15:28)] <- unname(c(
  EqAfterInvDonor, EqAfterInvDonorBulk))

MyData[CurrentIteration, c(29:51)] <- unname(c(
  EigValEq, DomEigVal, SignDomEigVal, SignEigValEqual,
  ComplexEigVal, EigValEqBulk, DomEigValBulk,
  SignDomEigValBulk, SignEigValEqualBulk,
  ComplexEigValBulk, smallchange, Mytmax, Mytstep))


BackupMyData <- MyData
MyData <- as.data.frame(MyData)
write.csv(MyData, file = paste0(DateTimeStamp, "outputdeSolveChangedContact", ".csv"),
          quote = FALSE, row.names = FALSE)


X <- limitseigenvalues[1]
ndigits <- 2
round(X - (10^-ndigits)/2, ndigits)

roundflex <- function(X, ndigits = 3, dir = "up") {
  if(dir == "down") {
    X <- round(X - (10^-ndigits)/2, ndigits)
  }
  if(dir == "up") {
    X <- round(X + (10^-ndigits)/2, ndigits)
  }
  return(X)
}

## An alternative way to highlight where equilibrium was not reached, but I could
# not figure out how to include info on the rectangles in the legend, and the
# creation of the raster dataframe requires re-entering the variable names, so
# is more error-prone. On the frames to indicate specific tiles, see
# https://stackoverflow.com/questions/13258454/marking-specific-tiles-in-geom-tile-geom-raster

frames = MyData[MyData$time==MyData$tmax, c("kp", "log10kn", "cd")]
ggplot(data = MyData, aes(x = kp, y = log10kn, fill = log10(PlasmidsEq))) + 
  geom_rect(colour = "white") + 
  scale_fill_gradient(low = "red", high = "green") +
  geom_rect(data=frames, size=1, fill=NA, colour="black",
            aes(xmin=kp - 0.05, xmax=kp + 0.05, ymin=log10kn - 0.05, ymax=log10kn + 0.05)) +
  # scale_colour_manual(name = "Ooh look", values = c("black"), labels = "Something cool") +
  facet_grid(cd ~ .) +
  labs(x = "kp",
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
ggplot(data = NULL, aes(x = kp, y = log10kn)) + 
  geom_tile(data = MyData, colour = "white") + 
  scale_fill_gradient(low = "red", high = "green") +
  geom_point(data = MyData[MyData$time==MyData$tmax, ], aes(colour = "black")) +
  scale_colour_manual(name = "", values = "black", labels = "Equilibrium not reached") +
  facet_grid(cd ~ .) +
  labs(x = "kp",
       y = "log10kn",
       title = "pair-formation model", 
       subtitle = "facet by cd ~ .") +
  theme(legend.position = "bottom")

# element_text(angle = 45, hjust = 1, vjust = 1, size = 14)
# theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1, vjust = 1))

#### Version history ####
# Before 9 January 2020 the equation for dR contained - kp*R*R instead of - 2*kp*R*R
# On 10 february 2020 I changed the equation of the nutrients to include washout of nutrients
# Later I changed kp, kn, gd, and gt to log10kp, log10kn, log10gd, and log10gt
# In version 9 (13 March 2020) I removed the biomass equations from the models. This has the
# advantage that the integration will be slightly faster since one equation less has to be
# integrated. More importantly, removing the biomass equation is needed to enable
# stability-analysis based on the eigenvalues of the jacobimatrix in the plasmid-free equilibrium
# On 13 March I also removed the plotting commands from the script and changed the structure of
# the output dataframe, to prevent repetition of the biomass and DInit = plasmid-bearing bacteria at the initial state.
# For some versions there is a script available using rootSolve to run fast to equilibrium,
# e.g. PairFormationSimulation2rootSolve8, but no plot over time can be made and it does not always work,
# and the output structure is different from deSolve, leading to different selection styles for storing data
# in the dataframe.
# In version 13 I corrected the calculation of the eigenvalues to use the plasmid-free equilibrium EqFull, not
# the plasmid-free equilibrium with DInit donors added to it (state). I also changed the calculation of the
# biomass at the initial state in the dataframe to a cleaner implementation.
# Version 13 is a copy of version 12, where I copied the script part by part to get the correct outlining
# of all the brackets, and I changed some comments and the order of some parts to get a clearer and cleaner
# script.
# In version 13 (?) added dots in the graphs to indicate the points where equilibrium was not reached.
# Besides, I splitted cSet and gSet in separate sets for D and T, and added DInit as a set and adjusted the
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
