#### Pair-formation model of conjugation ####


#### To do ####

# Also return ComplexEigVal and ComplexEigValBulk from function CalcEigenvalues

# Find method to obtain the order in which I specify the state from EstConjBulkDonor
# or from ModelEstConjBulkDonor, to prevent hardcoding names on the returned object DataEstConjBulk

# aes_string is soft-deprecated (see ?aes_string), use tidy evaluation idioms instead,
# see the quasiquotation section in aes() documentation.
# See also # On aes_string see https://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r

# See ggplot2::expand_limits to have the same limits and colorscale for the two plots

# Create function to filter on small negative and small positive state values and set them to 0

#  Make better comparison to check for unstable equilibrium where invasion is not possible,
#  currently comparison ofvalues that will not be exact

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
# last value is visible.

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


#### Plotting options ####
saveplots <- 0
tmaxsteady <- 1e8
timesParmsEst <- seq(from = 0, to = 3, by = 0.1)
MyColorBrew <- rev(brewer.pal(11, "Spectral")) # examples: display.brewer.all()


#### Functions ####
# Calculate the plasmid-free equilibrium (R*, Nutr*)
CalcEqPlasmidfree <- function(MyData) {
  with(as.list(MyData), {
    REq <- ((NI - w / bR)) / NutrConv
    NutrEq <- w / bR
    Eq <- c(NutrEq = NutrEq, REq = REq)
    print(Eq)
    print(class(Eq))
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
# adjusted pair-formation models for a short tail(timesParmsEst, 1) hours and
# calculate approximations of gdbulk and gtbulk from the output, following
# Zhong's approach for the calculations.
EstConjBulk <- function(MyData) {
  with(as.list(MyData), {
    state <- c(D = MyData[["DInit"]], R = MyData[["REq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesParmsEst, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)

        state <- c(R = MyData[["REq"]], Trans = MyData[["DInit"]], Mrt = 0, Mtt = 0)
    DataEstConjBulkTrans <- tail(ode(t = timesParmsEst, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)

    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
  })
}

# ODE-model describing pair-formation and conjugation. Pair-formation between
# plasmid-free recipients and plasmid-bearing donors or transconjugants depends
# on attachment rate kp. Conjugation from the donor or transconjugant occurs in
# the Mdr and Mrt pairs with intrinsic conjugation rates gd and gt respectively.
# This leads to formation of Mdt and Mtt pairs. The pairs fall apart with
# detachment rate kn. This structure of pair-formation is based on Zhong's model
# (Zhong 2010). I expanded the model to include costs, washout, and nutrients.
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
    dNutr <- (NI - Nutr)*w - NutrConv*Nutr*((1 - cd)*bR*D + bR*R + (1 - ct)*bR*Trans)
    dD <- (1 - cd)*bR*Nutr*D - w*D
    dR <- bR*Nutr*R - gdbulk*D*R - gtbulk*Trans*R - w*R
    dTrans <- (1 - ct)*bR*Nutr*Trans + gdbulk*D*R + gtbulk*Trans*R - w*Trans
    return(list(c(dNutr, dD, dR, dTrans)))
  })
}

# Numerically estimate the Jacobian matrix of the plasmid-free equilibrium of
# the models, then calculate (or approximate?) the eigenvalues of this matrix.
# The maximum real part of the eigenvalues is used to determine stability.
CalcEigenvalues <- function(MyData) {
  parms <- MyData
  EqFull <- c(Nutr = MyData[["NutrEq"]], D = 0, R = MyData[["REq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
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

# The initial state is the plasmid-free equilibrium (R*, Nutr*) with the
# addition of DInit donor bacteria per mL. Note that stol is based on the average
# of absolute rates of change, not the sum.
SimulationPairs <- function(InputSimulationPairs) {
  parms <- InputSimulationPairs
  state <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]],
             R = parms[["REq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelPairsNutr, parms = parms, stol = 1.25e-6)
  EqAfterInvDonor <- c(time = attr(out, "time"), out$y)
  return(EqAfterInvDonor)
}

# Run the bulk-conjugation model
SimulationBulk <- function(InputSimulationBulk) {
  parms <- InputSimulationBulk
  state <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]],
             R = parms[["REq"]], Trans = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelBulkNutr, parms = parms, stol = 2.5e-6)
  EqAfterInvDonor <- c(time = attr(out, "time"), out$y)
  return(EqAfterInvDonor)
}

# Create heatmaps, save if needed
CreatePlot <- function(fillvar, gradient2 = 0, limits = NULL, midpoint = 0, dataplot = MyData,
                       xvar = "log10(kp)", yvar = "log10(kn)", save = saveplots) {
  if(exists("DateTimeStamp") == FALSE) {
    warning("DateTimeStamp created to include in plot but does not correspond to filename of the dataset")
    DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")
  }
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar)) + 
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
    labs(caption = DateTimeStamp) +
    theme(legend.position = "bottom", plot.caption = element_text(vjust = 20))
  if(gradient2 == 1) {
    p <- p + scale_fill_gradient2(midpoint = midpoint)
  } else {
    p <- p + scale_fill_gradientn(colours = MyColorBrew, limits = limits)
    #p <- p + scale_fill_distiller(palette = "Spectral", direction = 1, limits = limits)
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
      print(paste("Taking the absolute value of the", nvaluesnegative, "values smaller than 0"), quote = FALSE)
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
# FileName <- "2020_juni_26_09_07_33outputsimulationcomplete.csv"
# MyData <- read.csv(FileName, header = TRUE, sep = ",", quote = "\"",
#                    dec = ".", stringsAsFactors = FALSE)
# MyData <- as.data.frame(MyData)
# DateTimeStamp <- substr(FileName, 1, nchar(FileName) - 28)

## Stable co-existence of recipients, transconjugants, and Mrt and Mtt pairs
# bRSet <- 1.7; NISet <- 10; kpSet <- 10^-9.6; knSet <- 10^0.5
# NutrConv <- 1E-6; wSet <- 0.04; cdSet <- 0.05; ctSet <- 0.05;
# gdSet <- 10^1.176; gtSet <- 10^1.176; DInitSet <- 1000
# Leads to time = 843129.4, Nutr = 0.02427729 D = 0, R = 3829217, Trans = 6142815,
# Mdr = Mdt = 0, Mrt = 324.6586, Mtt = 1520.435, timeBulk = 728451.6,
# NutrBulk = 0.02430818, DBulk = 0, RBulk = 3583808, TransBulk = 6391884.

## Small parameterset for tests
DInitSet <- c(1000)
bRSet <- c(1.7)
NISet <- c(10, 100)
NutrConv <- c(1E-6)
wSet <- c(0.04)
kpSet <- 10^seq(-10, -6, 2)
knSet <- 10^seq(-1, 3, 2)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- 15
gtSet <- 15

## Large dataset for tests
DInitSet <- c(500, 1E3)
bRSet <- c(0.8, 1.7)
NISet <- c(10, 100)
NutrConv <- 1e-6
wSet <- c(0.04, 0.06)
kpSet <- 10^seq(from = -11, to = -5, by = 0.5)
knSet <- 10^seq(from = -1, to = 3, by = 0.5)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(10, 15)
gtSet <- c(10, 15)

## This set worked (26 june 2020) (takes 8 minutes to run)
DInitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
NutrConv <- 1e-6
wSet <- c(0.04)
kpSet <- 10^seq(from = -10, to = -6, by = 0.25)
knSet <- 10^seq(from = -1, to = 3, by = 0.25)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(10, 15)
gtSet <- c(10, 15) 

## Try smaller steps for kp and kn (24 june 2020):
# works for pair-formation model, not for bulk-conjugation model
DInitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
NutrConv <- 1e-6
wSet <- c(0.04)
kpSet <- 10^seq(from = -10, to = -6, by = 0.2)
knSet <- 10^seq(from = -1, to = 3, by = 0.2)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(10, 15)
gtSet <- c(10, 15) 

## Using steps of 0.1 for kp and kn did not work
# Error: DLSODE- at T (=R1) and step size H (=R2), the corrector convergence failed repeatedly
# or with ABS(H)=HMIN. IN above message, R = 1.499002e+05, 2.091563e-09
# DLSODE- ISTATE illegal, in above message I=-5
# LSODE- run aborted, apparent infinite loop
# 
# Retry using jactype = "sparse", could also try using stode(s?) and/or supply jacobian


#### Main script ####

## Calculate plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(bR = bRSet, NI = NISet, NutrConv = NutrConv, w = wSet)
if(any(MyData <= 0)) warning("All parameters should have positive values.")

Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree))
Eqplasmidfree
MyData <- cbind(MyData, Eqplasmidfree)

if(any(Eqplasmidfree[, "REq"] <= 0)) {
  warning("The number of recipients at equilibrium is not positive!
Increase the nutrient concentration in the inflowing liquid by changing NI?")
}

print("Plasmid-free equilibrium calculated:")
print(Sys.time())

## Add combinations with the parameters needed to approximate gdbulk and gtbulk
MyData <- expand_grid(MyData, kp = kpSet, kn = knSet, gd = gdSet, gt = gtSet,
                      DInit = DInitSet)
DataEstConjBulk <- t(apply(X = MyData, MARGIN = 1, FUN = EstConjBulk)) # 97.47    0.00   97.51 

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulk <- MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor)
gdbulk <- unname(gdbulk)
MyData <- cbind(MyData, gdbulk = gdbulk)

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] + DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] + DataEstConjBulk[, "TransMrt"]
gtbulk <- MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans)
gtbulk <- unname(gtbulk)
MyData <- cbind(MyData, gtbulk = gtbulk)

print("Bulk-conjugation rates estimated:")
print(Sys.time())

MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)

MyInfoEigVal <- t(apply(MyData, MARGIN = 1, FUN = CalcEigenvalues))
MyData <- cbind(MyData, MyInfoEigVal)
print("Eigenvalues estimated:")
print(Sys.time())

DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")
write.csv(MyData, file = paste0(DateTimeStamp, "outputnosimulation.csv"),
          quote = FALSE, row.names = FALSE)


#### Plotting output ####

CreatePlot(fillvar = "NutrEq")
CreatePlot(fillvar = "REq")

# Create limits to have the same limits and colorscale for the two plots
limitsbulkrates <- c(floor(min(log10(c(MyData$gdbulk, MyData$gtbulk)))),
                     ceiling(max(log10(c(MyData$gdbulk, MyData$gtbulk)))))

CreatePlot(fillvar = "log10(gdbulk)", limits = limitsbulkrates)
CreatePlot(fillvar = "log10(gtbulk)", limits = limitsbulkrates)

PlotData <- MyData[which(MyData[, "cd"]==0.01 & MyData[, "ct"]==0.01 ), ]
CreatePlot(data = PlotData, fillvar = "pmax(gdbulk, gtbulk)/pmin(gdbulk, gtbulk)")
CreatePlot(data = PlotData, fillvar = "gdbulk/gtbulk")

CreatePlot(fillvar = "SignDomEigVal", gradient2 = 1)
CreatePlot(fillvar = "SignDomEigValBulk", gradient2 = 1)
CreatePlot(fillvar = "SignDomEigVal / SignDomEigValBulk", gradient2 = 1)

summary(MyData$SignDomEigVal - MyData$SignDomEigValBulk)
summary(MyData$SignDomEigVal / MyData$SignDomEigValBulk)

# Note: hardcoded legend
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
  labs(caption = DateTimeStamp) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("-1" = "darkblue", "1" = "darkred"),
                    name = "Stability",
                    labels = c("stable", "unstable"))
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "outputfactor(SignDomEigVal).png"))
}

# Note: hardcoded legend
ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal - SignDomEigValBulk))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
  labs(caption = DateTimeStamp) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("0" = "darkblue"),
                    name = "Difference in sign eigenvalues")
if(saveplots == 1 ) {
  ggsave(paste0(DateTimeStamp, "Difference in sign eigenvalues.png"))
}

# Create limits to have the same limits and colorscale for the two plots
DomEigVals <- c(MyData$DomEigVal, MyData$DomEigValBulk)
limitseigenvalues <- log10(range(DomEigVals[DomEigVals > 0]))
limitseigenvalues <- c(floor(limitseigenvalues[1]), ceiling(limitseigenvalues[2]))

CreatePlot(fillvar = "log10(DomEigVal)", limits = limitseigenvalues)
CreatePlot(fillvar = "log10(DomEigValBulk)", limits = limitseigenvalues)
CreatePlot(fillvar = "DomEigVal/DomEigValBulk")
CreatePlot(fillvar = "DomEigVal-DomEigValBulk")

summary(MyData$DomEigVal)
summary(MyData$DomEigValBulk)
summary(MyData$DomEigVal - MyData$DomEigValBulk)
summary(MyData$DomEigVal / MyData$DomEigValBulk)
summary(DomEigVals)

print("Finished plotting:")
print(Sys.time())

##################
BackupMyData <- MyData
# Now run pair-formation model for those cases where the dominant eigenvalue is negative

# If invasion is possible, run simulation to see how many bacteria of each
# population are present at equilibrium
IndexSimulation <- which(MyData$SignDomEigVal != -1)
ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))
InputSimulationPairs <- MyData[IndexSimulation, ColumnsToSelect]
OutputSimulationPairs <- t(apply(X = InputSimulationPairs, MARGIN = 1, FUN = SimulationPairs))

if(length(IndexSimulation) < nrow(MyData)) {
  NoSimulationNeeded <- cbind(time = 0, Nutr = MyData[-IndexSimulation, "NutrEq"],
                              D = 0, R = MyData[-IndexSimulation, "REq"],
                              Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  MyData <- rbind(cbind(MyData[IndexSimulation, ], OutputSimulationPairs),
                  cbind(MyData[-IndexSimulation, ], NoSimulationNeeded))
} else {
  MyData <- cbind(MyData[IndexSimulation, ], OutputSimulationPairs)
}

print("Pair-formation model completed running:")
print(Sys.time())
write.csv(MyData, file = paste0(DateTimeStamp, "outputsimulationpairs.csv"),
          quote = FALSE, row.names = FALSE)

IndexSimulationBulk <- which(MyData$SignDomEigValBulk != -1)
InputSimulationBulk <- MyData[IndexSimulationBulk, ColumnsToSelect]
OutputSimulationBulk <- t(apply(X = InputSimulationBulk, MARGIN = 1, FUN = SimulationBulk))
colnames(OutputSimulationBulk) <- paste0(colnames(OutputSimulationBulk), "Bulk")

if(length(IndexSimulationBulk) < nrow(MyData)) {
  NoSimulationNeededBulk <- cbind(timeBulk = 0, NutrBulk = MyData[-IndexSimulationBulk, "NutrEq"],
                                  DBulk = 0, RBulk = MyData[-IndexSimulationBulk, "REq"],
                                  TransBulk = 0)
  MyData <- rbind(cbind(MyData[IndexSimulationBulk, ], OutputSimulationBulk),
                  cbind(MyData[-IndexSimulationBulk, ], NoSimulationNeededBulk))
} else {
  MyData <- cbind(MyData[IndexSimulationBulk, ], OutputSimulationBulk)
}

print("Bulk-conjugation model completed running:")
print(Sys.time())

BackupMyData2 <- MyData

write.csv(MyData, file = paste0(DateTimeStamp, "outputsimulationpairsandbulk.csv"),
          quote = FALSE, row.names = FALSE)


MyData <- cbind(MyData, TotalD = NA, TotalR = NA, TotalTrans = NA, TotalPlasmid = NA, TotalBio = NA,
                TotalPlasmidBulk = NA, TotalBioBulk = NA)
MyData[, "TotalD"] <- MyData[, "D"] + MyData[, "Mdr"] + MyData[, "Mdt"]
MyData[, "TotalR"] <- MyData[, "R"] + MyData[, "Mdr"] + MyData[, "Mrt"]
MyData[, "TotalTrans"] <- MyData[, "Trans"] + MyData[, "Mdt"] + MyData[, "Mrt"] + 2*MyData[, "Mtt"]
MyData[, "TotalPlasmid"] <- MyData[, "TotalD"] + MyData[, "TotalTrans"]
MyData[, "TotalBio"] <- MyData[, "D"] + MyData[, "R"] + MyData[, "Trans"] + 2*MyData[, "Mdr"] +
  2*MyData[, "Mdt"] + 2*MyData[, "Mrt"] + 2*MyData[, "Mtt"]
MyData[, "TotalPlasmidBulk"] <- MyData[, "DBulk"] + MyData[, "TransBulk"]
MyData[, "TotalBioBulk"] <- MyData[, "DBulk"] + MyData[, "RBulk"] + MyData[, "TransBulk"]


write.csv(MyData, file = paste0(DateTimeStamp, "outputsimulationcomplete.csv"),
          quote = FALSE, row.names = FALSE)

names(MyData)
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
any(MyData$TotalR == MyData$TotalBio & MyData$SignDomEigVal == 1)
any(MyData$RBulk == MyData$TotalBioBulk & MyData$SignDomEigValBulk == 1)

# If this differs, comparisons between the 2 models should use fractions, not cell counts
summary(MyData$TotalBio / MyData$TotalBioBulk)
range(MyData$TotalBio / MyData$TotalBioBulk)

# If the biomass is the same accross simulations, one could use the number of cells instead of fractions of total biomass
summary(MyData$TotalBio)
max(MyData$TotalBio) / min(MyData$TotalBio) 
summary(MyData$TotalBioBulk)
max(MyData$TotalBioBulk) / min(MyData$TotalBioBulk)

# Which proportion of cells is plasmid-bearing, if it is not 0
min((MyData$TotalPlasmid / MyData$TotalBio)[which(MyData$TotalPlasmid / MyData$TotalBio != 0)])
max((MyData$TotalPlasmid / MyData$TotalBio)[which(MyData$TotalPlasmid / MyData$TotalBio != 0)])

min((MyData$TotalPlasmidBulk / MyData$TotalBioBulk)[which(MyData$TotalPlasmidBulk / MyData$TotalBioBulk != 0)])
max((MyData$TotalPlasmidBulk / MyData$TotalBioBulk)[which(MyData$TotalPlasmidBulk / MyData$TotalBioBulk != 0)])

A <- MyData[which(MyData[, "TotalPlasmid"]/MyData[, "TotalBio"] > 1e-3 & MyData[, "TotalPlasmid"]/MyData[, "TotalBio"] < 0.999), ]
dim(A)

RunAgain <- A[, c(1:15)]


# Compare biomass accross the two models
CreatePlot(fillvar = "TotalBio / TotalBioBulk", gradient2 = 1)

# Fraction plasmid-bearing cells
CreatePlot(fillvar = "TotalPlasmid/TotalBio")
CreatePlot(fillvar = "TotalPlasmidBulk/TotalBioBulk")
CreatePlot(fillvar = "(TotalPlasmid/TotalBio) / (TotalPlasmidBulk/TotalBioBulk)")

# Fraction donors
CreatePlot(fillvar = "TotalD/TotalBio")
CreatePlot(fillvar = "DBulk/TotalBioBulk")
CreatePlot(fillvar = "(DBulk/TotalBioBulk) / (TotalD/TotalBio)")

# Fraction recipients
CreatePlot(fillvar = "TotalR/TotalBio")
CreatePlot(fillvar = "RBulk/TotalBioBulk")
CreatePlot(fillvar = "(RBulk/TotalBioBulk) / (TotalR/TotalBio)")

# Fraction transconjugants
CreatePlot(fillvar = "TotalTrans/TotalBio")
CreatePlot(fillvar = "TransBulk/TotalBioBulk")
CreatePlot(fillvar = "(TransBulk/TotalBioBulk) / (TotalTrans/TotalBio)")

################################################################################

###### To create plots over time
# NOTE: CURRENTLY BROKEN (because no state is specified).

# Settings for simulations, plotting, and printing
mylty <- c(lty = c(3, 1, 2, 1, 1, 1, 1, 1))
# mycol <- c("black", "purple", "hotpink", "red", "yellow", "green1", "blue", "cyan")
mycol <- c("black", brewer.pal(7, "Set1"))
myylim <- c(1E-4, 1E7) # Defining the limits for the y-axis
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
plotoutput <- 0
extinctionthreshold <- 1E-10 # Population size is set to 0 if it is below the extinctionthreshold
verbose <- 0 # if verbose == 1, diagnositics on the simulations are printed and roots are indicated in the graphs
smallchange <- c(1E-5)
Mytmax <- c(1E5)
Mytstep <- c(10)
runsimulation <- 1 # if 0, no simulation is run, and only the
# parametervalues, plasmid-free equilibrium, and the eigenvalues are stored

# If runsimulation != 1 only the parametervalues, plasmid-free equilibrium, and the eigenvalues are stored

#### Create matrix to store data ####
Mydf <- expand.grid(DInit = DInitSet, bR = bRSet, NI = NISet, w = wSet, kpSet = kpSet,
                    knSet = knSet, cdSet = cdSet, ctSet = ctSet,
                    gdSet = gdSet, gtSet = gtSet, KEEP.OUT.ATTRS = FALSE)

TotalIterations <- length(bRSet)*length(NISet)*length(kpSet)*
  length(log10knSet)*length(wSet)*length(DInitSet)*length(cdSet)*length(ctSet)*
  length(gdSet)*length(gtSet)
print(TotalIterations)
CurrentIteration <- 0
MyData <- matrix(data = NA, nrow = TotalIterations, ncol = 61, byrow = TRUE) # To store output data
DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")

# Times for which output of the simulation is wanted.
# Note that I don't use timestep untill t = 100, and note furthermore that
# the used ode-solvers are variable-step methods, so the times in times
# are NOT the only times at which integration is performed. See
# ?diagnostics.deSolve() and ?lsodar() for details.
if(runsimulation == 1) {
  times <- c(0:100, seq(from = 100 + Mytstep, to = Mytmax, by = Mytstep)) 
}

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


out2 <- ode(t = times, y = state, func = ModelPairsNutr, parms = parmspair,
            rootfun = rootfun, events = list(func = eventfun, root = TRUE,
                                             terminalroot = 1))
if(verbose == TRUE) {
  print(diagnostics(out2))
  print(attributes(out2))
}
EqAfterInvDonor <- tail(out2, 1)[, ]

if(plotoutput == 1) {
  maintitle <- c("Pair-formation model")
  subtitlepair <- paste0("bR=", bRValue, " NI=", NIValue,
                         " log10kp=", log10kpValue, " log10kn=", log10knValue,
                         " e=", NutrConv, " w=", wValue,
                         " cd=", cdValue, " ct=", ctValue,
                         " log10gd=", log10gdValue, " log10gt=", log10gtValue)
  if(verbose == TRUE) {
    matplot.deSolve(out2, main = maintitle, sub = subtitlepair,
                    ylim = myylim, xlim = c(0, tail(attributes(out2)$troot, 2)[1]),
                    log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
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
} # End of preparing and plotting


if(runsimulation == 1) {
  # Run bulk-model
  out2bulk <- ode(t = times, y = stateBulk, func = ModelBulkNutr,
                  parms = parmsBulk, rootfun = rootfunBulk,
                  events = list(func = eventfunBulk, root = TRUE,
                                terminalroot = 1))
  if(verbose == TRUE) {
    print(diagnostics(out2bulk))
    print(attributes(out2bulk))
  }
  EqAfterInvDonorBulk <- tail(out2bulk, 1)[, ]
  
  if(plotoutput == 1) {
    subtitlebulk <- paste0("bR=", bRValue, " NI=", NIValue,
                           " log10kp=", log10kpValue, " log10kn=", log10knValue,
                           " e=", NutrConv, " w=", wValue,
                           " cd=", cdValue, " ct=", ctValue,
                           " log10gd=", log10gdValue, " log10gt=", log10gtValue,
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

# png(filename = paste0(DateTimeStamp, "longrun", CurrentIteration + 1, ".png"))
# dev.off()


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
