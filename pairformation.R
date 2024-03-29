#### Pair-formation model of conjugation ####

#### Introduction ####
# The pair-formation model was taken from equation (2) of Zhong (2010) and
# expanded to to include plasmid costs, washout, and nutrients.


#### References ####

# Zhong 2010: Zhong X, Krol JE, Top EM, Krone SM. 2010. Accounting for mating
# pair formation in plasmid population dynamics. Journal of Theoretical Biology
# 262:711-719.


#### To do ####

# aes_string is soft-deprecated (see help(aes_string)), use tidy evaluation idioms instead,
# see the quasiquotation section in aes() documentation and https://www.tidyverse.org/blog/2018/07/ggplot2-tidy-evaluation/
# See also # On aes_string see https://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r


#### Loading packages ####
# All packages are available from CRAN (https://cran.r-project.org/).
library(deSolve) # Solving differential equations with output over time.
library(dplyr) # filter(), near()
library(ggplot2) # Creating plots.
library(rootSolve) # Integrating ODEs, obtaining Jacobian matrix.
library(tidyr) # expand_grid() which allows for dataframe as input
library(cowplot) # plot_grid() to create a single figure from multiple plots,
                 # get_legend() to create shared legends


#### Plotting and simulation options ####
saveplots <- TRUE
plotdataapproxbulk <- FALSE
tmaxEstConj <- 3
tstepEstConj <- 0.1
timesEstConj <- seq(from = 0, to = tmaxEstConj, by = tstepEstConj)


#### Functions ####
# Calculate the plasmid-free equilibrium (R*, Nutr*) using the solution to
# dR/dt = R*(bR*Nutr / (Ks + Nutr) - w) == 0,
# dNutr/dt = (NI - Nutr)*w - NutrConv*Nutr*R*bR / (Ks + Nutr) == 0
# with conditions R > 0 and Nutr > 0
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
    if(plotdataapproxbulk == TRUE) {
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
    if(plotdataapproxbulk == TRUE) {
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

  EqFullBulk <- c(Nutr = MyData[["NutrEq"]], D = 0, R = MyData[["REq"]], Trans = 0)
  EigValBulk <- eigen(x = jacobian.full(y = EqFullBulk, func = ModelBulkNutr, parms = parms),
                      symmetric = FALSE, only.values = TRUE)$values
  ComplexEigValBulk <- is.complex(EigValBulk) 
  EigValBulk <- Re(EigValBulk)
  names(EigValBulk) <- paste0("EigvalBulk", 1:length(EigValBulk))
  DomEigValBulk <- max(EigValBulk)
  SignDomEigValBulk <- sign(DomEigValBulk)

  InfoEigVal <- c(EigVal, ComplexEigVal = ComplexEigVal, DomEigVal = DomEigVal,
                  SignDomEigVal = SignDomEigVal,
                  EigValBulk, ComplexEigValBulk = ComplexEigValBulk,
                  DomEigValBulk = DomEigValBulk,
                  SignDomEigValBulk = SignDomEigValBulk)
  return(InfoEigVal)
}

# Function to create and save heatmaps
CreatePlot <- function(dataplot = MyData, xvar = "log10(kp)", yvar = "log10(kn)",
                       fillvar = "factor(SignDomEigVal)",
                       filltype = "discrete", limits = NULL, 
                       labx = "Log10(attachment rate)",
                       laby = "Log10(detachment rate)",
                       filltitle, filllabels = c("No", "Yes"),
                       mytag = NULL,
                       manualvalues = c("TRUE" = "yellow", "FALSE" = "navy"),
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
    p <- p + scale_fill_viridis_d(filltitle,
                                  limits = if(is.null(limits)) {
                                    as.factor(c(-1, 1))
                                  } else {
                                    factor(limits)
                                  },
                                  labels = filllabels)
  }
  if(filltype == "continuous") {
    p <- p + scale_fill_viridis_c(filltitle, limits = limits)
  }
  if(filltype == "manual") {
    p <- p + scale_fill_manual(values = manualvalues, name = filltitle)
  }
  if(save == TRUE) {
    if(is.null(filename)) {
      fillvarname <- gsub("/", ".", fillvar)
      fillvarname <- gsub(" ", "", fillvarname)
      filename <- paste0(DateTimeStamp, fillvarname, ".png")      
    }
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      ggsave(filename)
    }
  }
  return(p)
}


#### Parameter values ####

## Explanation of parameters
# bR (1/h): maximum recipient growth rate
# w (1/h): washout rate
# Ks (mg/h): half-saturation constant for limiting nutrient 
# NI (mg/mL): nutrient concentration in the inflowing liquid  
# NutrConv (mg/cell): resource conversion rate
# DInit: the number of donors present at t=0 for simulations over time.
# kp (mL/(cell*h)): attachment rate
# kn (1/h): detachment rate
# gd (1/h): intrinsic conjugation rate from donor
# gt (1/h): intrinsic conjugation rate from transconjugant
# cd (unitless): costs in donor as fraction of recipient growth rate
# ct (unitless): costs in transconjugant as fraction of recipient growth rate
# gdbulk (mL/(cell*h)): bulk-conjugation rate of the donor
# gtbulk (mL/(cell*h)): bulk-conjugation rate of the transconjugant

## To read data from csv-file, uncomment this section and supply the 
# needed FileName
# FileName <- "YYYY_MM_DD_hh_mm.csv"
# MyData <- read.csv(FileName, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE)
# MyData <- as.data.frame(MyData)
# DateTimeStamp <- substr(FileName, 1, 16)


## Parameter set 1: show influence of nutrient concentration at the inflow and the
# washout rate on the stability of the plasmid-free equilibrium. The facet with
# the default values is also saved separately to compare the one-compartment
# pair-formation model and the two-compartment pair-formation model.
bRSet <- 0.738
wSet <- c(-log(0.5)/(24*10), -log(0.5)/24, -log(0.5)/(24/10))
Ks <- 0.004
NISet <- c(0.14, 1.4, 14)
NutrConvSet <- 1.4e-7
DInitSet <- 1000
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
knSet <- 10^seq(from = -1, to = 3, by = 0.1)
gdSet <- 15
gtSet <- 15
cdSet <- 0.18
ctSet <- 0.09

parameterset <- "1b"
NISet <- 10^seq(from = log10(0.14), to = log10(14), length.out = 9)
gtSet <- c(1, 15)
ctSet <- c(0.09, 0.18)

## Parameter set 2: show influence of costs, conjugation, attachment, and
# detachment rates on the stability of the plasmid-free equilibrium and on
# bulk-conjugation rates.
bRSet <- 0.738
wSet <- -log(0.5)/24
Ks <- 0.004
NISet <- 1.4
NutrConvSet <- 1.4e-7
DInitSet <- 1000
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
knSet <- 10^seq(from = -1, to = 3, by = 0.1)
gdSet <- c(1, 15)
gtSet <- c(1, 15)
cdSet <- c(0.09, 0.18)
ctSet <- c(0.09, 0.18)

## Parameter set 3: compare bulk- and pair-formation model over time
bRSet <- 0.738
wSet <- -log(0.5)/24
Ks <- 0.004
NISet <- 1.4
NutrConvSet <- 1.4e-7
DInitSet <- 1000
kpSet <- 10^c(-11, -9)
knSet <- 10
gdSet <- 15
gtSet <- 15
cdSet <- 0.18
ctSet <- 0.09


#### Run simulations ####

CheckParms <- c(bRSet, wSet, Ks, NISet, NutrConvSet, DInitSet, kpSet, knSet,
                gdSet, gtSet, cdSet, ctSet)
if(any(CheckParms <= 0)) {warning("All parameters should have positive values.")}
if(any(c(cdSet, ctSet) <= 0 | c(cdSet, ctSet) >= 1)) {
  warning("Costs should be larger than 0 and smaller than 1.")
}

TotalIterations <- length(bRSet)*length(wSet)*length(Ks)*length(NISet)*
  length(NutrConvSet)*length(DInitSet)*length(kpSet)*length(knSet)*
  length(gdSet)*length(gtSet)*length(cdSet)*length(ctSet)
print(paste(TotalIterations, "iterations to run."))

## Calculate plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(bR = bRSet, w = wSet, Ks = Ks, NI = NISet,
                      NutrConv = NutrConvSet)
Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree))
MyData <- cbind(MyData, Eqplasmidfree)

if(any(Eqplasmidfree <= 0)) {
  RowsNegativeEq <- sort(unique(which(Eqplasmidfree <= 0, arr.ind = TRUE)[, 1]))
  ColnamesNegativeEq <- colnames(Eqplasmidfree)[unique(which(Eqplasmidfree <= 0,
                                                             arr.ind = TRUE)[, 2])]
  warning("Plasmid-free equilibrium contains non-positive values in column(s) '",
          paste(ColnamesNegativeEq, collapse = "' and '"),
          "'.\nThe data will be included in the calculations anyway!
  This concerns the following rows of the dataframe: ",
          paste(RowsNegativeEq, collapse = ", "))
}

print(paste("Plasmid-free equilibrium calculated:", Sys.time()))

## Add combinations with the parameters needed to approximate gdbulk and gtbulk
MyData <- expand_grid(MyData, DInit = DInitSet, kp = kpSet, kn = knSet,
                      gd = gdSet, gt = gtSet)
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
print(paste("Bulk-conjugation rates estimated:", Sys.time()))

# Approximate eigenvalues
MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)

MyInfoEigVal <- t(apply(MyData, MARGIN = 1, FUN = CalcEigenvalues))
MyData <- cbind(MyData, MyInfoEigVal)
if(any(MyData[, "ComplexEigVal"] != 0)) {
  warning("Some eigenvalues of the pair-formation model have an imaginary part.")
}
if(any(MyData[, "ComplexEigValBulk"] != 0)) {
  warning("Some eigenvalues of the bulk-model have an imaginary part.")
}
print(paste("Eigenvalues estimated:", Sys.time()))
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
write.csv(MyData, file = paste0(DateTimeStamp, "pairformation.csv"),
          quote = FALSE, row.names = FALSE)

# Create facet labels and labeller 'function'
labw <- paste0("Washout rate ", signif(wSet, 2), "/h")
names(labw) <- wSet
labNI <- paste(signif(NISet, 2), "mg nutrients/mL\nat inflow")
names(labNI) <- NISet
labgd <- paste0("Donor:\nconjugation rate\n", gdSet, "/h")
names(labgd) <- gdSet
labgt <- paste0("Transconjugant:\nconjugation rate\n", gtSet, "/h")
names(labgt) <- gtSet
labcd <- paste("costs", cdSet)
names(labcd) <- cdSet
labct <- paste("costs", ctSet)
names(labct) <- ctSet
mylabeller <- labeller(w = labw, NI = labNI, gd = labgd, gt = labgt,
                       cd = labcd, ct = labct, .multi_line = FALSE,
                       .default = label_value)


#### Plotting output for parameter set 1 ####

# Influence of washout rate and nutrient concentration in the inflowing liquid
# on stability of the plasmid-free equilibrium (Figure 2 in article).
CreatePlot(filltitle = "Plasmid can invade", facetx = "NI", facety = "w")

CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk model)",
           facetx = "NI", facety = "w")

# Show if the largest eigenvalues of the pair-formation model and bulk model
# have equal signs
CreatePlot(fillvar = "factor(as.numeric(near(SignDomEigVal, SignDomEigValBulk)))",
           limits = c(0, 1),
           filltitle = "Possibility for invasion equal in pair-\nformation and bulk-conjugation models",
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
CreatePlot(fillvar = "log10(gdbulk)", filltype = "continuous",
           filltitle = "Log10(Donor bulkrate)",
           facetx = "NI", facety = "w")
CreatePlot(fillvar = "log10(gtbulk)", filltype = "continuous",
           filltitle = "Log10(Transconjugant\nbulkrate)",
           facetx = "NI", facety = "w")

# Create separate plot of the default facet to compare the one-compartment
# pair-formation model and the two-compartment pair-formation model (Figure 6A
# in article).
CreatePlot(dataplot = filter(MyData, near(w, -log(0.5)/24), near(NI, 1.4)),
           filltitle = "Plasmid can invade", facetx = ".", facety = ".", mytag = "A",
           filename = paste0(DateTimeStamp, "Figure6A.png"))

# Show percentage and counts of parameter combinations for which invasion is,
# or is not, possible
filteredDf <- NULL
for(i in NISet) {
  for(j in wSet) {
    MyDataFiltered <- filter(MyData, near(NI, i), near(w, j))
    invasion_n <- length(which(near(MyDataFiltered[, "SignDomEigVal"], 1)))
    no_invasion_n <- length(which(near(MyDataFiltered[, "SignDomEigVal"], -1)))
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
write.csv(filteredDf, file = paste0(DateTimeStamp, "invpercpar1.csv"),
          quote = FALSE, row.names = FALSE)


#### Plotting output for parameter set 1b ####
if(exists("parameterset") && parameterset == "1b") {
  filteredDf2 <- as_tibble(MyData) %>%
    group_by(gt, ct, NI, w, kn) %>%
    summarise(
      invasion_n = length(which(near(SignDomEigVal, 1))),
      no_invasion_n = length(which(near(SignDomEigVal, -1))),
      total_n = invasion_n + no_invasion_n,
      invasion_perc = round(100*invasion_n/total_n, 0),
      no_invasion_perc = round(100*no_invasion_n/total_n, 0),
      .groups = "drop"
    )
  write.csv(filteredDf2, file = paste0(DateTimeStamp, "invpercpar1Alternative1.csv"),
            quote = FALSE, row.names = FALSE)
  CreatePlot(dataplot = filteredDf2, xvar = "log10(NI)", yvar = "log10(kn)",
             fillvar = "invasion_perc", filltype = "continuous",
             labx = "Log10(nutrient concentration in inflow)", laby = "Log10(detachment rate)",
             filltitle = "Percentage of attachment rates for\nwhich invasion possible",
             facetx = "gt + ct", facety = "w", marginx = rep(1, 4),
             filename = paste0(DateTimeStamp, "Figure3.png"))
}


#### Plotting output for parameter set 2 ####

# To show influence of costs and intrinsic conjugation rates on stability of the
# plasmid-free equilibrium (Figure A5 in the supplement):
CreatePlot(filltitle = "Plasmid can invade")

# Stability of the equilibrium for the bulk-conjugation model
# (plot not shown in article)
CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk model)")

# Show if sign of largest eigenvalues for pair-formation and bulk model is the 
# same
CreatePlot(fillvar = "factor(as.numeric(near(SignDomEigVal, SignDomEigValBulk)))",
           limits = c(0, 1),
           filltitle = "Possibility for invasion equal in pair-\nformation and bulk-conjugation models")

## Compare simulation-derived and calculated bulk-rates
# Bulk-rates calculated as gdbulk = gd * kp / kn and gtbulk = gt * kp / kn.
# First filter data on conjugation rates and costs. Copy data on
# bulk-conjugation rates and eigenvalues from simulation-derived
# bulk-conjugation rates to new columns, and replace the old columns with data
# from calculated bulk-conjugation rates
DataCalcRates <- filter(MyData, near(gd, 15), near(gt, 15),
                        near(cd, 0.18), near(ct, 0.09))
DataCalcRates[, "gdbulksim"] <- DataCalcRates[, "gdbulk"]
DataCalcRates[, "gtbulksim"] <- DataCalcRates[, "gtbulk"]
DataCalcRates[, "SignDomEigValsim"] <- DataCalcRates[, "SignDomEigVal"]
DataCalcRates[, "SignDomEigValBulksim"] <- DataCalcRates[, "SignDomEigValBulk"]
DataCalcRates[, "gdbulk"] <- DataCalcRates[, "gd"] * DataCalcRates[, "kp"] /
  DataCalcRates[, "kn"]
DataCalcRates[, "gtbulk"] <- DataCalcRates[, "gt"] * DataCalcRates[, "kp"] /
  DataCalcRates[, "kn"]

# Remove existing columns with data on eigenvalues, calculate eigenvalues using
# calculated bulk-conjugation rates, and save the results
DataCalcRates <- select(DataCalcRates, -c("SignDomEigVal", "SignDomEigValBulk",
                                          starts_with(c("Eigval", "ComplexEigVal", "DomEigVal"))))
MyInfoEigValCalcRates <- t(apply(DataCalcRates, MARGIN = 1, FUN = CalcEigenvalues))
DataCalcRates <- cbind(DataCalcRates, MyInfoEigValCalcRates)
DataCalcRates <- select(DataCalcRates, -starts_with("Eigval"))
write.csv(DataCalcRates, file = paste0(DateTimeStamp, "calcbulkratepairformation.csv"),
          quote = FALSE, row.names = FALSE)

## Plots to compare the influence of conjugation, attachment, and detachment
# rates on the simulation-derived and calculated bulk-conjugation rates
# Plots of data from simulation-derived bulk-rates
limitscalcbulkrates <- range(log10(c(DataCalcRates$gtbulk, DataCalcRates$gtbulksim)))
theme_no_y <- theme(legend.position = "none", axis.title.y = element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0), "pt"))

Plotgtbulksim <- CreatePlot(dataplot = DataCalcRates,
                            fillvar = "log10(gtbulksim)", filltype = "continuous",
                            limits = limitscalcbulkrates,
                            filltitle = "Log10(Transconjugant\nbulkrate)",
                            facetx = ".", facety = ".", save = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
PlotLegendgtbulk <- get_legend(Plotgtbulksim)
Plotgtbulksim <- Plotgtbulksim + theme(legend.position = "none")
Plotgtbulkcalc <- CreatePlot(dataplot = DataCalcRates,
                             fillvar = "log10(gtbulk)", filltype = "continuous",
                             limits = limitscalcbulkrates, filltitle = NULL,
                             facetx = ".", facety = ".", save = FALSE) + theme_no_y
Plotgtbulk <- plot_grid(Plotgtbulksim, Plotgtbulkcalc) +
  theme(plot.margin = unit(c(10, 0, 10, 0), "pt"))

PlotSignEigValsim <- CreatePlot(dataplot = DataCalcRates,
                                fillvar = "factor(SignDomEigValBulksim)",
                                filltitle = "Plasmid can invade\n(bulk model)",
                                facetx = ".", facety = ".", save = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
PlotLegendSignEigVal <- get_legend(PlotSignEigValsim)
PlotSignEigValsim <- PlotSignEigValsim + theme(legend.position = "none")
PlotSignEigValcalc <- CreatePlot(dataplot = DataCalcRates,
                                 fillvar = "factor(SignDomEigValBulk)",
                                 filltitle = NULL,
                                 facetx = ".", facety = ".", save = FALSE) + theme_no_y
PlotSignEigVal <- plot_grid(PlotSignEigValsim, PlotSignEigValcalc) +
  theme(plot.margin = unit(c(10, 0, 10, 0), "pt"))

PlotSignEigValEqualsim <- CreatePlot(dataplot = DataCalcRates,
                                     fillvar = "factor(as.numeric(near(SignDomEigValsim, SignDomEigValBulksim)))",
                                     limits = c(0, 1),
                                     filltitle = "Possibility for invasion equal in\nbulk- and pair-formation models",
                                     facetx = ".", facety = ".", save = FALSE) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
PlotLegendSignEigValEqual <- get_legend(PlotSignEigValEqualsim)
PlotLegendSignEigValEqual <- PlotLegendSignEigValEqual
PlotSignEigValEqualsim <- PlotSignEigValEqualsim + theme(legend.position = "none")
PlotSignEigValEqualcalc <- CreatePlot(dataplot = DataCalcRates,
                                      fillvar = "factor(as.numeric(near(SignDomEigVal, SignDomEigValBulk)))",
                                      limits = c(0, 1), filltitle = NULL,
                                      facetx = ".", facety = ".", save = FALSE) + theme_no_y
PlotSignEigValEqual <- plot_grid(PlotSignEigValEqualsim, PlotSignEigValEqualcalc) +
  theme(plot.margin = unit(c(10, 0, 10, 0), "pt"))

# Collect all sub-plots in one figure (Figure A4 in the supplement)
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
png(filename = paste0("FigA4", DateTimeStamp, ".tiff"),
    width = 420, height = 1.5*480)
plot_grid(Plotgtbulk, PlotLegendgtbulk,
          PlotSignEigVal, PlotLegendSignEigVal,
          PlotSignEigValEqual, PlotLegendSignEigValEqual,
          ncol = 1, byrow = TRUE,
          rel_heights = c(1, 0.1, 1, 0.1, 1, 0.1),
          labels = c("A", "", "B", "", "C"), hjust = -1) +
  theme(plot.margin = margin(0, 0, 5, 0, "pt"))
dev.off()


# show plots comparing bulk-conjugations derived from simulations with calculated
# bulk-conjugation rates.
CreatePlot(dataplot = DataCalcRates, fillvar = "gdbulk / gdbulksim",
           filltype = "continuous", filltitle = "gdbulkcalc / gdbulksim",
           facetx = ".", facety = ".")
CreatePlot(dataplot = DataCalcRates, fillvar = "log10(gdbulk / gdbulksim)",
           filltype = "continuous", filltitle = "Log10(gdbulkcalc / gdbulksim)",
           facetx = ".", facety = ".")
CreatePlot(dataplot = DataCalcRates, fillvar = "gtbulk / gtbulksim",
           filltype = "continuous", filltitle = "gtbulkcalc / gtbulksim",
           facetx = ".", facety = ".")
CreatePlot(dataplot = DataCalcRates, fillvar = "log10(gtbulk / gtbulksim)",
           filltype = "continuous", filltitle = "Log10(gtbulkcalc / gtbulksim)",
           facetx = ".", facety = ".")


# Change layout of labels for next plots
labgd <- paste0("Donor conjugation rate ", gdSet, "/h")
names(labgd) <- gdSet
labgt <- paste0("Transconjugant conjugation rate ", gtSet, "/h")
names(labgt) <- gtSet
mylabeller <- labeller(gd = labgd, gt = labgt, .default = label_value)

# The influence of conjugation, attachment, and detachment rates on the
# bulk-conjugation rates (Figure 7 in the article and Figure S1 in the supplement).
# Data is filtered to contain only one value for costs, because costs do not
# influence the bulk-conjugation rates.
limitsbulkrates <- range(log10(c(MyData$gdbulk, MyData$gtbulk)))
CreatePlot(dataplot = filter(MyData, near(gd, 15) &
                               near(cd, cdSet[1]) & near(ct, ctSet[1])),
           fillvar = "log10(gtbulk)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant bulkrate)",
           facetx = "gt", facety = ".")

CreatePlot(dataplot = filter(MyData, near(gt, 15) &
                               near(cd, cdSet[1]) & near(ct, ctSet[1])),
           fillvar = "log10(gdbulk)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate)",
           facetx = "gd", facety = ".")

filteredDf <- NULL
for(k in cdSet) {
  for(l in gdSet) {
    for(i in ctSet) {
      for(j in gtSet) {
        MyDataFiltered <- filter(MyData, near(ct, i), near(gt, j),
                                 near(cd, k), near(gd, l))
        invasion_n <- length(which(near(MyDataFiltered[, "SignDomEigVal"], 1)))
        no_invasion_n <- length(which(near(MyDataFiltered[, "SignDomEigVal"], -1)))
        total_n <- invasion_n + no_invasion_n
        invasion_perc <- round(100*invasion_n/total_n, 0)
        no_invasion_perc <- round(100*no_invasion_n/total_n, 0)
        filteredDf <- rbind(filteredDf,
                            data.frame(cd = k, gd = l, ct = i, gt = j,
                                       invasion_perc = invasion_perc,
                                       no_invasion_perc = no_invasion_perc,
                                       invasion_n = invasion_n,
                                       no_invasion_n = no_invasion_n,
                                       total_n = total_n))
      }
    }
  }
}
print(filteredDf)
write.csv(filteredDf, file = paste0(DateTimeStamp, "invpercpar2.csv"),
          quote = FALSE, row.names = FALSE)


#### Plotting output for parameter set 3 ####

# Note:
# Generating the file names with the DateTimeStamp as used above led to problems
# because multiple plots were generated within one second. They got the same
# name and were not all saved. I now use a more precise DateTimeStamp. For
# alternative to get sequential numbers, see remark in ?ggsave() on use of
# filename = "figure%03d.png".

# Settings for simulations, plotting, and printing
mylty <- c(lty = c(3, 1, 2, 1, 1, 1, 1, 1))
mycol <- c("black", "blue", "red", "darkgreen", "brown", "purple", "darkorange",
           "yellow")
myylim <- c(1E-6, 1E7) # Defining the limits for the y-axis
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
plotoutput <- 1
verbose <- 0 # if verbose == 1, diagnostics on the simulations are printed
Mytmax <- c(4000)
Mytstep <- c(10)

## Select subset of data to use for runs over time
ColumnsToSelect <- c(1:(which(names(MyData) == "Eigval1") - 1))
MyDataSubset <- MyData[, ColumnsToSelect]
print(paste(dim(MyDataSubset)[1], "iterations to run."))

# Times for which output of the simulation is wanted. The used ODE-solvers are
# variable-step methods, so the times in times are NOT the only times at which
# integration is performed. See help(diagnostics.deSolve()) and help(lsodar())
# (both in the deSolve package) for details.
times <- c(0:100, seq(from = 100 + Mytstep, to = Mytmax, by = Mytstep))

# Note: the functions RunOverTime and PlotOverTime in the two-compartment script are
# more elaborate to enable comparison of D, R, Trans in the bulk-conjugation model
# with TotalD, TotalR, TotalTrans in the pair-formation model. That might be a
# fairer comparison.
RunOverTime <- function(parms = MyDataSubset, verbose = FALSE, printsubtitle = FALSE, ...) {
  state <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]], R = parms[["REq"]],
             Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  out2 <- ode(t = times, y = state, func = ModelPairsNutr, parms = parms,
              verbose = verbose)
  EqAfterInvasion <- tail(out2, 1)
  print(EqAfterInvasion)
  if(verbose == TRUE) {
    print(diagnostics(out2))
    print(attributes(out2))
  }
  # Change variable name to be same as in equations of the manuscript
  colnames(out2)[which(colnames(out2) == "Mdr")] <- "Mrd"
  PlotOverTime(plotdata = out2, parms = parms, type = "Pair", verbose = verbose,
               printsubtitle = printsubtitle, saveplot = saveplots)
  stateBulk <- c(Nutr = parms[["NutrEq"]], D = parms[["DInit"]],
                 R = parms[["REq"]], Trans = 0)
  out2bulk <- ode(t = times, y = stateBulk, func = ModelBulkNutr, parms = parms,
                  verbose = verbose)
  EqAfterInvasionBulk <- tail(out2bulk, 1)
  if(verbose == TRUE) {
    print(diagnostics(out2bulk))
    print(attributes(out2bulk))
  }
  PlotOverTime(plotdata = out2bulk, parms = parms, type = "Bulk",
               verbose = verbose, printsubtitle = printsubtitle,
               saveplot = saveplots)
  
  EqAfterInvasionTotal <- cbind(EqAfterInvasion, EqAfterInvasionBulk)
  names(EqAfterInvasionTotal) <- c(colnames(EqAfterInvasion),
                                   paste0(colnames(EqAfterInvasionBulk), "Bulk"))
  return(EqAfterInvasionTotal)
}

PlotOverTime <- function(plotdata = out2, parms = parms, type = "Pair",
                         printsubtitle = FALSE, saveplot = saveplots, ...) {
  maintitle <- paste(type, "model")
  subtitle <- paste0("bR=", parms[["bR"]], " NI=", parms[["NI"]],
                     " log10kp=", log10(parms[["kp"]]), " log10kn=", log10(parms[["kn"]]),
                     " NutrConv=", parms[["NutrConv"]], "\nw=", signif(parms[["w"]], digits = 4),
                     " cd=", parms[["cd"]], " ct=", parms[["ct"]])
  if(type == "Pair") {
    subtitle <- paste0(subtitle, " gd=", parms[["gd"]], " gt=", parms[["gt"]]) 
  } else {
    subtitle <- paste0(subtitle, " gdbulk=",  signif(parms[["gdbulk"]], digits = 4),
                       " gtbulk=", signif(parms[["gtbulk"]], digits = 4))
  }
  if(printsubtitle == FALSE) {
    subtitle <- NULL
  } 
  filename <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M_%OS3"),
                     "runovertime", gsub(" ", "", maintitle), ".png")
  if(saveplot == TRUE & file.exists(filename) == FALSE) {
    png(filename = filename)
    matplot.deSolve(plotdata, main = maintitle, sub = subtitle, ylab = "Density",
                    ylim = myylim, log = if(yaxislog == 1) {"y"}, col = mycol,
                    lty = mylty, lwd = 2, legend = list(x = "topright"))
    grid()
    dev.off()
  } else {
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    }
    matplot.deSolve(plotdata, main = maintitle, sub = subtitle, ylab = "Density",
                    ylim = myylim, log = if(yaxislog == 1) {"y"}, col = mycol,
                    lty = mylty, lwd = 2, legend = list(x = "topright"))
    grid()
  }
}

# Plots are shown in Figure A2 in the supplement
EqAfterInvasion <- t(apply(X = MyDataSubset, MARGIN = 1, FUN = RunOverTime,
                           printsubtitle = FALSE))
EqAfterInvasion <- cbind(MyDataSubset, EqAfterInvasion)
write.csv(EqAfterInvasion, file = paste0(DateTimeStamp, "runovertimeonecomp.csv"),
          quote = FALSE, row.names = FALSE)
