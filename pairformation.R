#### Pair-formation model of conjugation ####

#### To do ####

# Naming of objects inside EstConjBulk() should be updated, see DataEstConjBulk
# and EstConjBulk() in older versions of the script. Consider using mutate(),
# see the script of the two-compartment model.

# Use diff() in EstConjBulk() instead of manually calculating the difference
# between consecutive values

# If I want to include stochastics, see Price 'An efficient moments-based inference
# method for within-host bacterial infection dynamics' and Volkova 'Modelling dynamics
# of plasmid-gene mediated antimicrobial resistance in enteric bacteria using stochastic
# differential equations' for ideas

# Also return ComplexEigVal and ComplexEigValBulk from function CalcEigenvalues

# Instead of 'biomass', I should use 'density'.

# The calculations of TotalDEstConjBulkDonor, TotalREstConjBulkDonor, TotalTransEstConjBulkDonor
# ect. for approximating bulk-rates can be moved inside the EstConjBulk function.

# Terminology: residence time or retention time?

# Prevent overlapping values for the colorbar (could use angle with hjust, vjust)

## Voor plots over time waar ik het pair-formation en het bulk-model met elkaar Vergelijk
# is plotten van TotalD, TotalR, TotalTrans van het pair-model en DBulk, RBulk, TransBulk van het bulk-model
# een betere vergelijking, zie PlotOverTime uit het script voor het twee-compartimenten model.

# For PlotOverTime: use \n to print subtitle over two lines, use oma (see ?par and ?mtext)
# to prevent text from overlapping with legend

# Using CreatePlot(gradient2 = TRUE) geeft probleem dat factor(limits) gebruikt
# wordt terwijl default limits = NULL, dus
# scale_fill_viridis_d(limits = factor(limits)) veranderen in iets als 
# scale_fill_viridis_d(!if(is.null(limits)) {limits = factor(limits)})
# om dat te voorkomen?

# aes_string is soft-deprecated (see help(aes_string)), use tidy evaluation idioms instead,
# see the quasiquotation section in aes() documentation and https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
# See also # On aes_string see https://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r

# See ggplot2::expand_limits to have the same limits and colorscale for the two plots

# In functions use more informative name arguments instead of e.g. mydf = mydf.

# Make better comparison to check for unstable equilibrium where invasion is not
# possible, currently comparison of values that will not be exact. For example,
# use dplyr::near(SignDomEigVal, 1)

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

# The functions RunOverTime and PlotOverTime in the two-compartment script are
# more elaborate to enable comparison of D, R, Trans in the bulk-conjugation model
# with TotalD, TotalR, TotalTrans in the pair-formation model.
# In the two-compartment model I have not used root- and eventfunctions,
# maybe should also delete them here, since it sometimes takes very long to
# execute code (if populations go extinct and therefor many roots are found?)

# Try to get different plot objects together as different panels in one figure.
# See ggpubr:ggarrange(labels = c("A", "B"), common.legend = TRUE, legend = "bottom")

# Use vectors for atol, to have different tolerances for the cell-densities (~1),
# and nutrient concentration (~1*10^-8 ?)

# To compare the pair-formation and bulk-model select rows where any of the two
# models have instable equilibrium, i.e. which(pmax(MyData[, "SignDomEigVal"],
# MyData[, "SignDomEigValBulk"]) == 1)

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
library(RColorBrewer) # For better color schemes [NOTE: only used for plots over time]
library(rootSolve) # Integration, obtaining jacobian matrix and eigenvalues.
library(tidyr) # for 'expand.grid()' with dataframe as input
library(dplyr) # for mutate() to create new variables in dataframe/tibble

#### Plotting and simulation options ####
saveplots <- 1
atol <- 1e-10 # lower absolute error tolerance of integrator used by runsteady()
# to prevent 'DLSODE-  Warning..internal T (=R1) and H (=R2) are [1] 0 such that
# in the machine, T + H = T on the next step  [1] 0 (H = step size). Solver will
# continue anyway', which eventually leads to aborted integration.
tmaxsteady <- 1e8
tmaxEstConj <- 3
tstepEstConj <- 0.1
timesEstConj <- seq(from = 0, to = tmaxEstConj, by = tstepEstConj)
plotdataapproxbulk <- 0

#### Functions ####
# Calculate the plasmid-free equilibrium (R*, Nutr*) using the solution to
# dR/dt = R*(bR*Nutr / (Ks + Nutr) - w) == 0,
# dNutr/dt = (NI - Nutr)*w - NutrConv*Nutr*R*bR / (Ks + Nutr) ==
# 0 with R > 0 and Nutr > 0
CalcEqPlasmidfree <- function(MyData) {
  with(as.list(MyData), {
    NutrEq <- w*Ks / (bR - w)
    REq <- (NI - NutrEq) / NutrConv
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
    state <- c(D = MyData[["DInit"]], R = MyData[["REq"]], Trans = 0,
               Mdr = 0, Mdt = 0, Mrt = 0)
    parms <- MyData
    DataEstConjBulkDonor <- ode(t = timesEstConj, y = state,
                                func = ModelEstConjBulkDonor, parms = parms)
    if(plotdataapproxbulk == 1) {
      subtitle <- paste0("log10(kp,kn)=", log10(MyData[["kp"]]),
                         " ", log10(MyData[["kn"]]))
      matplot.deSolve(DataEstConjBulkDonor, ylim = c(1E-7, 1E7), log = "y",
                      col = c("black", "purple", "green1", "red", "yellow", "hotpink"),
                      lty = c(1, 2, 1, 1, 1, 1), lwd = 2,
                      legend = list(x = "bottomright"), sub = subtitle)
      grid()
    }
    DataEstConjBulkDonor <- tail(DataEstConjBulkDonor, 1)
    state <- c(R = MyData[["REq"]], Trans = MyData[["DInit"]], Mrt = 0, Mtt = 0)
    DataEstConjBulkTrans <- ode(t = timesEstConj, y = state,
                                func = ModelEstConjBulkTrans, parms = parms)
    if(plotdataapproxbulk == 1) {
      matplot.deSolve(DataEstConjBulkTrans, ylim = c(1E-7, 1E7), log = "y",
                      col = c("purple", "green1", "hotpink", "cyan"),
                      lty = c(2, 1, 1, 1), lwd = 2,
                      legend = list(x = "bottomright"), sub = subtitle)
      grid()
    }
    DataEstConjBulkTrans <- tail(DataEstConjBulkTrans, 1)
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
    dNutr <- (NI - Nutr)*w - NutrConv*Nutr*((1 - cd)*bR*(D + Mdr + Mdt) +
                              bR*(R + Mdr + Mrt) + (1 - ct)*bR*(Trans + Mdt + Mrt + 2*Mtt))/(Ks + Nutr)
    dD <- (1 - cd)*bR*Nutr*(D + Mdr + Mdt)/(Ks + Nutr) - kp*D*R + kn*(Mdr + Mdt) - w*D
    dR <- bR*Nutr*(R + Mdr + Mrt)/(Ks + Nutr) - kp*R*(D + Trans) + kn*(Mdr + Mrt) - w*R
    dTrans <- (1 - ct)*bR*Nutr*(Trans + Mdt + Mrt + 2*Mtt)/(Ks + Nutr) - kp*R*Trans +
      kn*(Mdt + Mrt + 2*Mtt) - w*Trans
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
    dNutr <- (NI - Nutr)*w - NutrConv*bR*Nutr*((1 - cd)*D + R + (1 - ct)*Trans)/(Ks + Nutr)
    dD <- ((1 - cd)*bR*Nutr/(Ks + Nutr) - w)*D
    dR <- (bR*Nutr/(Ks + Nutr) - gdbulk*D - gtbulk*Trans - w)*R
    dTrans <- (1 - ct)*bR*Nutr*Trans/(Ks + Nutr) + gdbulk*D*R + gtbulk*Trans*R - w*Trans
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

# Function to create and save heatmaps
CreatePlot <- function(dataplot = MyData, xvar = "log10(kp)", yvar = "log10(kn)",
                       fillvar = "factor(SignDomEigVal)",
                       filltype = "discrete", limits = NULL, 
                       labx = "Log10(attachment rate)",
                       laby = "Log10(detachment rate)",
                       filltitle, filllabels = c("No", "Yes"), mytag = NULL,
                       manualvalues = c("TRUE" = "darkgreen", "FALSE" = "red"),
                       facetx = "gt + ct", facety = "gd + cd", as.table = TRUE,
                       marginx = NULL, marginy = NULL,
                       save = saveplots, filename = NULL, addstamp = FALSE, ...) {
  if(addstamp == TRUE & exists("DateTimeStamp") == FALSE) {
    warning("DateTimeStamp created to include in plot does not correspond to filename of the dataset")
    DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
  }
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar),
              subtitle = subtitle) + 
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_fixed(ratio = 1, expand = FALSE) +
    facet_grid(as.formula(paste(facety, "~", facetx)), as.table = as.table,
               labeller = mylabeller) +
    theme(legend.position = "bottom") +
    labs(x = labx, y = laby, tag = mytag)
  if(!is.null(marginx)) {
    p <- p + theme(strip.text.x = element_text(margin = margin(marginx)))
  }
  if(!is.null(marginy)) {
    p <- p + theme(strip.text.y = element_text(margin = margin(marginy)))
  }
  if(addstamp == TRUE) {
    p <- p + labs(caption = DateTimeStamp) +
      theme(plot.caption = element_text(vjust = 20))
  }
  if(filltype == "discrete") {
    p <- p + scale_fill_viridis_d(filltitle, limits = if(is.null(limits)) {NULL
      } else {factor(limits)}, labels = filllabels)
  }
  if(filltype == "continuous") {
    p <- p + scale_fill_viridis_c(filltitle, limits = limits)
  }
  if(filltype == "manual") {
    p <- p + scale_fill_manual(values = manualvalues, name = filltitle)
  }
  print(p)
  if(save == TRUE) {
    if(is.null(filename)) {
      fillvarname <- gsub("/", ".", fillvar)
      fillvarname <- gsub(" ", "", fillvarname)
      filename <- paste0(DateTimeStamp, "output", fillvarname, ".png")      
    }
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      ggsave(filename)
    }
  }
}


#### Parameter values ####

## Explanation of parameters
# Growthrates (1/h) for recipient: bR (growth rates for donor and transconjugant
# are not explicitly defined, but calculated as (1 - cd)*bR and (1 - ct)*bR)
# Conjugation rates (1/h) from donor and transconjugant: gd and gt
# Mating pair attachment rate (mL * cell^-1 * h^-1): kp (k+ in Zhong's notation)
# Mating pair detachment rate (1/h): kn (k- in Zhong's notation)
# Nutrient concentration in the inflowing liquid (milligram * mL^-1): NI
# Resource conversion rate (milligram per cell division): NutrConv
# Washout rate (1/h): w

## To read data from csv-file
# FileName <- "2021_01_25_11_38outputnosimtwocomp.csv"
# MyData <- read.csv(FileName, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE)
# MyData <- as.data.frame(MyData)
# DateTimeStamp <- substr(FileName, 1, 16)


## Parameterset 1: show NI and w influence stability of plasmid-free equilibrium
DInitSet <- c(1E3)
bRSet <- c(0.738)
Ks <- 0.004
NISet <- c(0.14, 1.4, 14)
NutrConvSet <- 1.4e-7
# Washout rates corresponding to median residence times of 240, 24, and 2.4 h
# (corresponding to mean residence times of 346, 35, and 3.5 hours).
wSet <- c(-log(0.5)/(24*10), -log(0.5)/24, -log(0.5)/(24/10))
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
knSet <- 10^seq(from = -1, to = 3, by = 0.1)
cdSet <- c(0.18)
ctSet <- c(0.09)
gdSet <- c(15)
gtSet <- c(15)

## Parameterset 2: show influence of costs, conjugation, attachment, and 
## detachment rates on the stability of the plasmid-free equilibrium and on
## bulk-conjugation rates.
DInitSet <- c(1E3)
bRSet <- c(0.738)
Ks <- 0.004
NISet <- 1.4
NutrConvSet <- 1.4e-7
wSet <- -log(0.5)/24
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
knSet <- 10^seq(from = -1, to = 3, by = 0.1)
cdSet <- c(0.09, 0.18)
ctSet <- c(0.09, 0.18)
gdSet <- c(1, 15)
gtSet <- c(1, 15)

## Parameterset 3: to compare bulk- and pair model over time
DInitSet <- c(1E3)
bRSet <- c(0.738)
Ks <- 0.004
NISet <- 1.4
NutrConvSet <- 1.4e-7
wSet <- -log(0.5)/24
kpSet <- 10^c(-12, -10, -8)
knSet <- 10
cdSet <- c(0.18)
ctSet <- c(0.09)
gdSet <- 15
gtSet <- 15


#### Main script ####

CheckParms <- c(DInitSet, bRSet, NISet, Ks, NutrConvSet, wSet, kpSet, knSet, cdSet, ctSet)
if(any(CheckParms <= 0)) warning("All parameters should have positive values.")
if(any(c(cdSet, ctSet) >= 1)) warning("Costs should be larger than 0 and smaller than 1.")

TotalIterations <- length(DInitSet)*length(bRSet)*length(NISet)*length(Ks)*length(NutrConvSet)*
  length(wSet)*length(kpSet)*length(knSet)*length(cdSet)*length(ctSet)*
  length(gdSet)*length(gtSet)
TotalIterations

## Calculate plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(bR = bRSet, NI = NISet, Ks = Ks, NutrConv = NutrConvSet, w = wSet)

Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree))
MyData <- cbind(MyData, Eqplasmidfree)

if(any(Eqplasmidfree[, "REq"] <= 0)) {
  warning("Some rows contain non-positive recipient densities at equilibrium.
  They have been discarded. Increase growth rate or the nutrient concentration 
  in the inflowing liquid, or decrease the outflow rate to prevent this.")
  MyData <- MyData[-which(MyData[, "REq"] < 0), ]
}

print("Plasmid-free equilibrium calculated:")
print(Sys.time())

## Add combinations with the parameters needed to approximate gdbulk and gtbulk
MyData <- expand_grid(MyData, gd = gdSet, gt = gtSet, DInit = DInitSet,
                      kp = kpSet, kn = knSet)
DataEstConjBulk <- t(apply(X = MyData, MARGIN = 1, FUN = EstConjBulk))

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] +
  DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] +
  DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulk <- unname(MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] /
                   (TotalDEstConjBulkDonor * TotalREstConjBulkDonor))

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] +
  DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] +
  DataEstConjBulk[, "TransMrt"]
gtbulk <- unname(MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] /
                   (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans))

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

#### Create facet labels and labeller 'function' ####
labNI <- paste(signif(NISet, 2), "mg nutrients/mL\nat inflow")
names(labNI) <- NISet
labw <- paste0("Washout rate ", signif(wSet, 2), "/h")
names(labw) <- wSet
labcd <- paste("Donor costs", cdSet)
names(labcd) <- cdSet
labct <- paste("Transconjugant\ncosts", ctSet)
names(labct) <- ctSet
labgd <- paste0("Donor conjugation\nrate ", gdSet, "/h")
names(labgd) <- gdSet
labgt <- paste0("Transconjugant\nconjugation rate\n", gtSet, "/h")
names(labgt) <- gtSet

mylabeller <- labeller(NI = labNI, w = labw, cd = labcd, ct = labct,
                       gd = labgd, gt = labgt, .default = label_both)

#### Plotting output for parameterset 1 ####

# Influence of washout rate and nutrient concentration in the inflowing liquid
# on stability of the plasmid-free equilibrium (Figure 2 in article).
CreatePlot(filltitle = "Plasmid can invade", facetx = "NI", facety = "w")

# Show if the dominant eigenvalues of the pair-formation model and bulk model
# have equal signs (Figure S1 in article)
CreatePlot(fillvar = "factor(SignDomEigVal == SignDomEigValBulk)",
           filltype = "manual",
           filltitle = "Dominant eigenvalues\nhave equal signs",
           facetx = "NI", facety = "w")

CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk model)",
           facetx = "NI", facety = "w")

# Nutrient concentration at the inflow influences recipient cell density,
# whereas the washout rate influences nutrient concentration (not shown in
# article).
CreatePlot(fillvar = "log10(REq)", filltype = "continuous",
           filltitle = "Log10(Recipient density)",
           facetx = "NI", facety = "w")
CreatePlot(fillvar = "log10(NutrEq)", filltype = "continuous",
           filltitle = "Log10(Nutrient concentration)",
           facetx = "NI", facety = "w")

# Bulk conjugation rates
limitsbulkrates <- range(log10(c(MyData$gdbulk, MyData$gtbulk)))
CreatePlot(dataplot = filter(MyData, cd == cdSet[1] & ct == ctSet[1]),
           fillvar = "log10(gdbulk)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate)",
           facetx = "NI", facety = "w", save = FALSE)
CreatePlot(dataplot = filter(MyData, cd == cdSet[1] & ct == ctSet[1]),
           fillvar = "log10(gtbulk)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant bulkrate)",
           facetx = "NI", facety = "w", save = FALSE)

# Show percentage and counts of parameter combinations for which invasion is,
# or is not, possible
filteredDf <- NULL
for(i in NISet) {
  for(j in wSet) {
    MyDataFiltered <- filter(MyData, NI == i, w == j)
    invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == 1))
    no_invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == -1))
    total_n <- invasion_n + no_invasion_n
    invasion_perc <- round(100*invasion_n/total_n, 0)
    no_invasion_perc <- round(100*no_invasion_n/total_n, 0)
    
    filteredDf <- rbind(filteredDf,
                        data.frame(NI = i, w = j, invasion_perc = invasion_perc,
                                   no_invasion_perc = no_invasion_perc,
                                   invasion_n = invasion_n,
                                   no_invasion_n = no_invasion_n, total_n = total_n))
  }
}
print(filteredDf)

#### Plotting output for parameterset 2 ####

# To show influence of costs and intrinsic conjugation rates on stability of the
# plasmid-free equilibrium (Figure 3 in article):
CreatePlot(filltitle = "Plasmid can invade",
           marginx = c(0,0,0,0), marginy = c(0,0,0,0))

# Stability of the equilibrium for the bulk-conjugation model
# (plot not shown in article)
CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk model)",
           marginx = c(0,0,0,0), marginy = c(0,0,0,0))

# Show if sign of dominant eigenvalues for pair-formation and bulk model is the 
# same (Figure S2 in article)
CreatePlot(fillvar = "factor(SignDomEigVal == SignDomEigValBulk)",
           filltype = "manual",
           filltitle = "Dominant eigenvalues\nhave equal signs",
           marginx = c(0,0,0,0), marginy = c(0,0,0,0))

# Change layout of labels for next plots
labgd <- paste0("Donor conjugation rate ", gdSet, "/h")
names(labgd) <- gdSet
labgt <- paste0("Transconjugant conjugation rate ", gtSet, "/h")
names(labgt) <- gtSet
mylabeller <- labeller(NI = labNI, w = labw, cd = labcd, ct = labct,
                       gd = labgd, gt = labgt, .default = label_both)

## The influence of conjugation, attachment, and detachment rates on the 
## bulk-conjugation rates (Figure 4 in the article). Data is filtered to show
## only one value for costs, because costs do not influence the bulk-conjugation
## rates.
limitsbulkrates <- range(log10(c(MyData$gdbulk, MyData$gtbulk)))
CreatePlot(dataplot = filter(MyData, gt == 15 & cd == cdSet[1] & ct == ctSet[1]),
           fillvar = "log10(gdbulk)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate)",
           facetx = "gd", facety = "gt")

CreatePlot(dataplot = filter(MyData, gd == 15 & cd == cdSet[1] & ct == ctSet[1]),
           fillvar = "log10(gtbulk)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant bulkrate)",
           facetx = "gt", facety = "gd")

filteredDf <- NULL
for(k in cdSet) {
  for(l in gdSet) {
    for(i in ctSet) {
      for(j in gtSet) {
        MyDataFiltered <- filter(MyData, ct == i, gt == j, cd == k, gd == l)
        invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == 1))
        no_invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == -1))
        total_n <- invasion_n + no_invasion_n
        invasion_perc <- round(100*invasion_n/total_n, 0)
        no_invasion_perc <- round(100*no_invasion_n/total_n, 0)
        filteredDf <- rbind(filteredDf,
                            data.frame(ct = i, gt = j, invasion_perc = invasion_perc,
                                       no_invasion_perc = no_invasion_perc,
                                       invasion_n = invasion_n,
                                       no_invasion_n = no_invasion_n, total_n = total_n))
      }
    }
  }
}
print(filteredDf)


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
Mytmax <- c(1E4)
Mytstep <- c(10)

#### Create matrix to store data ####
TheseRows <- c(1, nrow(MyData))
TheseRows <- 1:nrow(MyData)
# TheseRows <- ceiling(seq(from = nrow(MyData)/(4*4*2*2), to = nrow(MyData), length.out = 4*4*2))
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
    # filename <- paste0(DateTimeStamp, "output", gsub(" ", "", maintitle), ".png")
    # See remark in ?ggsave() on use of filename = "figure%03d.png" to get sequential numbers.
    # filename <- paste0(DateTimeStamp, "output", gsub(" ", "", maintitle), "%03d", ".png")
    filename <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M_%OS3"), "output", gsub(" ", "", maintitle), ".png") 
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
    dev.off()
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

# Plots are shown in Figure S3 in the article.
EqAfterInvasion <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime))

EqAfterInvasion <- cbind(Mydf, EqAfterInvasion)

write.csv(EqAfterInvasion, file = paste0(DateTimeStamp, "outputrunovertimeonecomp.csv"),
          quote = FALSE, row.names = FALSE)
