#### Two-compartment pair-formation model of conjugation ####

#### To do ####
# Also see the 'To do' section in the pair-formation script for the one-compartment model

# Instead of checking if the plasmid-free equilibrium that has been found is valid,
# I could derive the invasion-criterion of the cell-free equilibrium
# (i.e. Nutr* = NI, D*, R*, ... = 0) to determine if cells can exist.

# Can I use the migration rates to set better states for the plasmid-free equilibrium

# Rethink the way of defining donor growth rate and costs: unless it is assumed 
# that the donor is the same species as the recipient, it does not make sense
# to define donor growth rates as recipient growth rate altered by costs from
# plasmid carriage. Instead, the donor could be assumed to have a growth rate
# unrelated to the recipient growth rate, and if the donor is adapted to the
# plasmid it does not have plasmid costs, but the new environment could be assumed
# to lead to costs in growth as well.

# See Macken 1994 'The dynamics of bacteria-plasmid systems' for analytic treatment
# and a way to get rid of the nutrient equation (but is that feasible if I want
# to modify nutrient levels over time?)

# Add CreatePlot() and SummaryPlot() function (see script for one-compartment
# model)

# Check if using MigrLumWall = MigrWallLum = 0 en DInitWall = 0 leads to same 
# results as the single-compartment model.

# Also use ModelEstConjBulkDonor to approximate bulk conjugation rates from the transconjugant
# (note: than Mdr should again be added in calculation of TotalREstConjBulkTrans)

# The calculations of TotalDEstConjBulkDonor, TotalREstConjBulkDonor, TotalTransEstConjBulkDonor
# ect. for approximating bulk-rates can be moved inside the EstConjBulk function ?

# In calculation of bulkrates use knLum, knWall instead of kn and knWall, then
# use kn = knLum as function argument for lumen and kn = knWall as function argument for wall.
# This will prevent creating 'DataWall' as separate dataframe

## !! !! POSSIBLE ERROR !! !!
# ! Should ScaleAreaPerVolume be the inverse (i.e. exchange multiplication and
# division by ScaleAreaPerVolume). This factor includes two aspects: (1) what is
# the ratio of volume and surface occupied by the bacteria (i.e., how much
# square cm of surface do the the bacteria contained in Y cubic cm occupy, see
# Imran (2005) for some discussion of this) and (2) what is the ratio of the
# volume and surface of the system that is modelled (i.e., if the gut is assumed
# to be an open square tube, the surface in square cm equals 4 times the
# volume in cubic cm). When modelling a round tube, the volume / surface = r / 2,
# irrespective of chosen length.

# Saving plots in PlotOverTime() (called from RunOverTime()) fails after the
# first plot because names are the same. Add 'Iteration' as column to data
# and append the value from that column to the name? That would als make it
# easier to see which plot is which iteration (so also print it in graph subtitles).

## To check if other parameters than kp, kn, gd, gt, cd, ct which are separately
# shown in the plot have multiple values (which leads to plotting multiple values
# on top of each other) the following can be used:
# # Parms <- list(DInitLumSet = DInitLumSet, DInitWallSet = DInitWallSet, bRSet = bRSet,
# NISet = NISet, NutrConvSet = NutrConvSet, MigrLumWallSet = MigrLumWallSet,
# MigrWallLumSet = MigrWallLumSet, ScaleAreaPerVolSet = ScaleAreaPerVolSet,
# wSet = wSet, kpSet = kpSet, knSet = knSet, kpWallSet = kpWallSet,
# knWallSet = knWallSet, cdSet = cdSet, ctSet = ctSet)
# ColIndexRelParms <- which(lengths(Parms) > 1)
# NamesRelParms <- names(ColIndexRelParms)
# expand.grid(Parms[ColIndexRelParms])
# NamesToSkip <- c("kpSet", "knSet", "cdSet", "ctSet", "gdSet", "gtSet")
# ColIndexRelParmsSelection <- ColIndexRelParms[-c(which(names(ColIndexRelParms) %in% NamesToSkip))]
# NamesRelParmsSelection <- names(ColIndexRelParmsSelection)
# ParmsSetsToSubset <- expand.grid(Parms[NamesRelParmsSelection])

#### References ####

# Zhong 2010: Zhong X, Krol JE, Top EM, Krone SM. 2010. Accounting for mating
# pair formation in plasmid population dynamics. Journal of Theoretical Biology
# 262:711-719.

# Imran 2005: Imran M, Jones D, Smith H. 2005. Biofilms and the plasmid
# maintenance question. Mathematical biosciences 193:183-204.

#### Loading packages ####
library(deSolve) # Solve differential equations with results over time.
library(ggplot2) # For plotting data
library(RColorBrewer) # For better color schemes
library(rootSolve) # Integration, obtaining jacobian matrix and eigenvalues.
library(tidyr) # for 'expand.grid()' with dataframe as input
library(dplyr) # for mutate() to create new variables in dataframe/tibble

#### Plotting options ####
saveplots <- 0
atol <- 1e-10 # lower absolute error tolerance of integrator used by runsteady()
# to prevent 'DLSODE-  Warning..internal T (=R1) and H (=R2) are [1] 0 such that
# in the machine, T + H = T on the next step  [1] 0 (H = step size). Solver will
# continue anyway', which eventually leads to aborted integration.
tmaxsteady <- 1e8
timesEstConj <- seq(from = 0, to = 3, by = 0.1)
MyColorBrew <- rev(brewer.pal(11, "Spectral")) # examples: display.brewer.all()

#### Parameter values ####
# MigrLumWallSet: rate of migration from the lumen to the wall
# MigrWallLumSet: rate of migration from the wall to the lumen
# See script of the one-compartment model for explanation of the other parameter values

# FileName <- "2020_september_09_12_18_07MyDataNoCarryCap.csv"
# MyData <- read.csv(FileName, header = TRUE, sep = ",", quote = "\"",
#                   dec = ".", stringsAsFactors = FALSE)
# MyData <- as.data.frame(MyData)
# DateTimeStamp <- substr(FileName, 1, nchar(FileName) - 36)

# Parametervalues
NILumSet <- c(1, 10, 100)
NIWallSet <- c(1, 10, 100)
wLumSet <- c(0.04)
wWallSet <- c(0.01)
NutrConvSet <- c(1e-6)
bRSet <- c(1.7)
MigrLumWallSet <- 0.03
MigrWallLumSet <- 0.05
ScaleAreaPerVolSet <- c(1)
DInitLumSet <- 1E3
DInitWallSet <- 1E3
cdSet <- 0.05
ctSet <- 0.01
kpSet <- 10^-8
kpWallSet <- 10^-8
# kpWallSet <- 10^seq(from = -10, to = -6, by = 1)
knSet <- 1
knWallSet <- 1
# knWallSet <- 10^seq(from = -1, to = 3, by = 1)
gdSet <- 15
gtSet <- 15

#### Functions ####
ModelBulkNutrPlasmidfree <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutrLum <- (NILum - NutrLum)*wLum - NutrConv*bR*NutrLum*RLum
    dRLum <- (bR*NutrLum - wLum - MigrLumWall)*RLum + MigrWallLum*RWall*ScaleAreaPerVol
    dNutrWall <- (NIWall - NutrWall)*wWall - NutrConv*bR*NutrWall*RWall
    dRWall <- (bR*NutrWall - wWall - MigrWallLum)*RWall + MigrLumWall*RLum/ScaleAreaPerVol
    return(list(c(dNutrLum, dRLum, dNutrWall, dRWall)))
  })
}

# I used the equilibrium values for the one-compartment model as initial state values.
RunToPlamidfreeEq <- function(parms) {
  with(as.list(parms), {
    state <- c(NutrLum = NutrLumGuess, RLum = RLumGuess,
               NutrWall = NutrWallGuess, RWall = RWallGuess)
    out <- runsteady(y = state, time = c(0, tmaxsteady),
                     func = ModelBulkNutrPlasmidfree, parms = parms,
                     stol = 1.25e-6, atol = atol)
    Eqplasmidfree <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
    names(Eqplasmidfree) <- paste0(names(Eqplasmidfree), "Init")
    return(Eqplasmidfree)
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

# Functions to estimate bulk-conjugation rates by running simulations with the
# adjusted pair-formation models for a short time (i.e., tail(timesEstConj, 1)
# hours) and calculate approximations of gdbulk and gtbulk from the output,
# following Zhong's approach for the calculations.
EstConjBulkLum <- function(MyData) {
  with(as.list(MyData), {
    if(MyData[["DInitLum"]] == 0) {
      warning("DInitLum == 0, DInitWall will be used to approximate bulkrates in the lumen instead!")
      state <- c(D = MyData[["DInitWall"]], R = MyData[["RLumInit"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    } else {
      state <- c(D = MyData[["DInitLum"]], R = MyData[["RLumInit"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    }
    
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    if(MyData[["DInitLum"]] == 0) {
      state <- c(R = MyData[["RLumInit"]], Trans = MyData[["DInitWall"]], Mrt = 0, Mtt = 0)
    } else {
      state <- c(R = MyData[["RLumInit"]], Trans = MyData[["DInitLum"]], Mrt = 0, Mtt = 0)
    }
    
    DataEstConjBulkTrans <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)
    
    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
  })
}

EstConjBulkWall <- function(MyData) {
  with(as.list(MyData), {
    if(MyData[["DInitWall"]] == 0) {
      warning("DInitWall == 0, DInitLum will be used to approximate bulkrates at the wall instead!")
      state <- c(D = MyData[["DInitLum"]], R = MyData[["RWallInit"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    } else {
      state <- c(D = MyData[["DInitWall"]], R = MyData[["RWallInit"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    }
    
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    if(MyData[["DInitWall"]] == 0) {
      state <- c(R = MyData[["RWallInit"]], Trans = MyData[["DInitLum"]], Mrt = 0, Mtt = 0)
    } else {
      state <- c(R = MyData[["RWallInit"]], Trans = MyData[["DInitWall"]], Mrt = 0, Mtt = 0)
    }
    
    DataEstConjBulkTrans <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)
    
    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
  })
}

# Bulk-conjugation model, with inflow, outflow, conversion by bacteria for nutrient
# equations, and growth, washout, migration from and to compartment, conjugation
# from donors and transconjugants for bacterial equations.
ModelBulkNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutrLum <- (NILum - NutrLum)*wLum - NutrConv*bR*NutrLum*((1 - cd)*DLum + RLum + (1 - ct)*TransLum)
    dDLum <- ((1 - cd)*bR*NutrLum - wLum - MigrLumWall)*DLum + MigrWallLum*DWall*ScaleAreaPerVol
    dRLum <- (bR*NutrLum - wLum - MigrLumWall)*RLum + MigrWallLum*RWall*ScaleAreaPerVol -
      (gdbulkLum*DLum + gtbulkLum*TransLum)*RLum
    dTransLum <- ((1 - ct)*bR*NutrLum - wLum - MigrLumWall)*TransLum + MigrWallLum*TransWall*ScaleAreaPerVol +
      (gdbulkLum*DLum + gtbulkLum*TransLum)*RLum
    dNutrWall <- (NIWall - NutrWall)*wWall - NutrConv*bR*NutrWall*((1 - cd)*DWall + RWall + (1 - ct)*TransWall)
    dDWall <- ((1 - cd)*bR*NutrWall - wWall - MigrWallLum)*DWall + MigrLumWall*DLum/ScaleAreaPerVol
    dRWall <- (bR*NutrWall - wWall - MigrWallLum)*RWall + MigrLumWall*RLum/ScaleAreaPerVol -
      (gdbulkWall*DWall - gtbulkWall*TransWall)*RWall
    dTransWall <- ((1 - ct)*bR*NutrWall - wWalll - MigrWallLum)*TransWall + MigrLumWall*TransLum/ScaleAreaPerVol +
      (gdbulkWall*DWall + gtbulkWall*TransWall)*RWall
    return(list(c(dNutrLum, dDLum, dRLum, dTransLum, dNutrWall, dDWall, dRWall, dTransWall)))
  })
}

# ODE-model describing pair-formation and conjugation for the two-compartment
# model, including migration between the compartments. Pair-formation between
# plasmid-free recipients and plasmid-bearing donors or transconjugants depends
# on attachment rates kp in the lumen and kpWall at the wall. Conjugation from
# donors or transconjugants occurs in the Mdr and Mrt pairs with intrinsic
# conjugation rates gd and gt respectively. This leads to formation of Mdt and
# Mtt pairs. The pairs fall apart with detachment rates kn in the lumen and
# knWall at the wall. The structure of pair-formation is based on Zhong's model
# (Zhong 2010). I expanded the model to include costs in growth for
# plasmid-bearing bacteria, washout from the lumen, and nutrients.
# The structure of the two compartments and migration between them is based on
# Imran (2005), but I added a carrying capacity in the wall-compartment limiting
# migration from the lumen to the wall.
ModelPairsNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    BioWall <- DWall + RWall + TransWall + MdrWall + MdtWall + MrtWall + MttWall
    
    dNutr <- (NI - Nutr)*w -
      NutrConv*Nutr*(
        (1 - cd)*bR*(DLum + MdrLum + MdtLum + ScaleAreaPerVol*(DWall + MdrWall + MdtWall)) +
          bR*(RLum + MdrLum + MrtLum + ScaleAreaPerVol*(RWall + MdrWall + MrtWall)) +
          (1 - ct)*bR*(TransLum + MdtLum + MrtLum + 2*MttLum +
                         ScaleAreaPerVol*(TransWall + MdtWall + MrtWall + 2*MttWall))
      )
    
    dDLum <- (1 - cd)*bR*Nutr*(DLum + MdrLum + MdtLum) - kp*DLum*RLum +
      kn*(MdrLum + MdtLum) - (w + MigrLumWall*(1 - BioWall/KWall))*DLum + MigrWallLum*DWall*ScaleAreaPerVol
    dRLum <- bR*Nutr*(RLum + MdrLum + MrtLum) - kp*RLum*(DLum + TransLum) +
      kn*(MdrLum + MrtLum) - (w + MigrLumWall*(1 - BioWall/KWall))*RLum + MigrWallLum*RWall*ScaleAreaPerVol
    dTransLum <- (1 - ct)*bR*Nutr*(TransLum + MdtLum + MrtLum + 2*MttLum) -
      kp*RLum*TransLum + kn*(MdtLum + MrtLum + 2*MttLum) -
      (w + MigrLumWall*(1 - BioWall/KWall))*TransLum + MigrWallLum*TransWall*ScaleAreaPerVol
    dMdrLum <- kp*DLum*RLum - kn*MdrLum - gd*MdrLum - (w + MigrLumWall*(1 - BioWall/KWall))*MdrLum +
      MigrWallLum*MdrWall*ScaleAreaPerVol
    dMdtLum <- gd*MdrLum - kn*MdtLum - (w + MigrLumWall*(1 - BioWall/KWall))*MdtLum +
      MigrWallLum*MdtWall*ScaleAreaPerVol
    dMrtLum <- kp*RLum*TransLum - kn*MrtLum - gt*MrtLum -
      (w + MigrLumWall*(1 - BioWall/KWall))*MrtLum + MigrWallLum*MrtWall*ScaleAreaPerVol
    dMttLum <- gt*MrtLum - kn*MttLum -
      (w + MigrLumWall*(1 - BioWall/KWall))*MttLum + MigrWallLum*MttWall*ScaleAreaPerVol
    
    dDWall <- (1 - cd)*bR*Nutr*(DWall + MdrWall + MdtWall) - kpWall*DWall*RWall +
      knWall*(MdrWall + MdtWall) + MigrLumWall*(1 - BioWall/KWall)*DLum/ScaleAreaPerVol - MigrWallLum*DWall
    dRWall <- bR*Nutr*(RWall + MdrWall + MrtWall) - kpWall*RWall*(DWall + TransWall) +
      knWall*(MdrWall + MrtWall) + MigrLumWall*(1 - BioWall/KWall)*RLum/ScaleAreaPerVol - MigrWallLum*RWall
    dTransWall <- (1 - ct)*bR*Nutr*(TransWall + MdtWall + MrtWall + 2*MttWall) -
      kpWall*RWall*TransWall + knWall*(MdtWall + MrtWall + 2*MttWall) +
      MigrLumWall*(1 - BioWall/KWall)*TransLum/ScaleAreaPerVol - MigrWallLum*TransWall
    dMdrWall <- kpWall*DWall*RWall - knWall*MdrWall - gd*MdrWall +
      MigrLumWall*(1 - BioWall/KWall)*MdrLum/ScaleAreaPerVol - MigrWallLum*MdrWall
    dMdtWall <- gd*MdrWall - knWall*MdtWall + MigrLumWall*(1 - BioWall/KWall)*MdtLum/ScaleAreaPerVol - 
      MigrWallLum*MdtWall
    dMrtWall <- kpWall*RWall*TransWall - knWall*MrtWall - gt*MrtWall + 
      MigrLumWall*(1 - BioWall/KWall)*MrtLum/ScaleAreaPerVol - MigrWallLum*MrtWall
    dMttWall <- gt*MrtWall - knWall*MttWall + 
      MigrLumWall*(1 - BioWall/KWall)*MttLum/ScaleAreaPerVol - MigrWallLum*MttWall
    
    return(list(c(dNutr, dDLum, dRLum, dTransLum, dMdrLum, dMdtLum, dMrtLum, dMttLum,
                  dDWall, dRWall, dTransWall, dMdrWall, dMdtWall, dMrtWall, dMttWall)))
  })
}

# Numerically estimate the Jacobian matrix of the plasmid-free equilibrium of
# the models, then calculate (or approximate?) the eigenvalues of this matrix.
# The maximum real part of the eigenvalues is used to determine stability.
CalcEigenvalues <- function(MyData) {
  parms <- MyData
  EqFull <- c(Nutr = MyData[["NutrInit"]], DLum = 0, RLum = MyData[["RLumInit"]],
              TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0, MttLum = 0,
              DWall = 0, RWall = MyData[["RWallInit"]],
              TransWall = 0, MdrWall = 0, MdtWall = 0, MrtWall = 0, MttWall = 0)
  EigVal <- eigen(x = jacobian.full(y = EqFull, func = ModelPairsNutr, parms = parms),
                  symmetric = FALSE, only.values = TRUE)$values
  ComplexEigVal <- is.complex(EigVal) 
  EigVal <- Re(EigVal)
  names(EigVal) <- paste0("Eigval", 1:length(EigVal))
  DomEigVal <- max(EigVal)
  SignDomEigVal <- sign(DomEigVal)
  SignEigValEqual <- identical(rep(SignDomEigVal, length(EigVal)), sign(Re(EigVal)))
  
  EqFullBulk <- c(Nutr = MyData[["NutrInit"]], DLum = 0, RLum = MyData[["RLumInit"]], TransLum = 0, 
                  DWall = 0, RWall = MyData[["RWallInit"]], TransWall = 0)
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
# equilibrium (Nutr*, RLum*, RWall*) with the addition of DInitLum and DInitWall
# donor bacteria per mL to the lumen and wall compartment, respectively is used
# as initial state. Note that stol is based on the average of absolute rates of
# change, not the sum.
SimulationPairs <- function(InputSimulationPairs) {
  parms <- InputSimulationPairs
  state <- c(Nutr = parms[["NutrInit"]], DLum = parms[["DInitLum"]],
             RLum = parms[["RLumInit"]], TransLum = 0, MdrLum = 0, MdtLum = 0,
             MrtLum = 0, MttLum = 0, DWall = parms[["DInitWall"]],
             RWall = parms[["RWallInit"]], TransWall = 0, MdrWall = 0, MdtWall = 0,
             MrtWall = 0, MttWall = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelPairsNutr,
                   parms = parms, stol = 1.25e-6, atol = atol)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

SimulationBulk <- function(InputSimulationBulk, state) {
  parms <- InputSimulationBulk
  state <- c(Nutr = parms[["NutrInit"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumInit"]], TransLum = 0, 
             DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]], TransWall = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelBulkNutr,
                   parms = parms, stol = 2.5e-6, atol = atol)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

CreatePlot <- function(fillvar, gradient2 = 0, limits = NULL, midpoint = 0,
                       dataplot = MyData, xvar = "log10(kp)", yvar = "log10(kn)",
                       facetx = "cd + gd", facety = "ct + gt", save = saveplots, ...) {
  if(exists("DateTimeStamp") == FALSE) {
    warning("DateTimeStamp created to include in plot but does not correspond to filename of the dataset")
    DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")
  }
  mycaption <- paste(DateTimeStamp)
  
  p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar), subtitle = subtitle) + 
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(as.formula(paste(facetx, "~", facety)), labeller = label_both) +
    labs(caption = mycaption) +
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

CreatePlot2 <- function(fillvar, gradient2 = 0, limits = NULL, midpoint = 0,
                        dataplot = MyData, xvar = "log10(kp)", yvar = "log10(kn)",
                        facetx = "cd + gd", facety = "ct + gt", save = saveplots, ...) {
  CumRowIndex <- NULL
  iteration <- 1
  dataplottotal <- dataplot
  for(kpWallsubset in unique(dataplottotal[, "kpWall"])) {
    for(knWallsubset in unique(dataplottotal[, "knWall"])) {
      subtitle <- paste0("kpWall= ", kpWallsubset, " knWall=", knWallsubset)
      RowIndex <- dataplottotal[, "kpWall"] == kpWallsubset & dataplottotal[, "knWall"] == knWallsubset
      dataplot <- dataplottotal[RowIndex, ]
      if(exists("DateTimeStamp") == FALSE) {
        warning("DateTimeStamp created to include in plot but does not correspond to filename of the dataset")
        DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")
      }
      mycaption <- paste(DateTimeStamp, subtitle)
      
      p <- ggplot(data = dataplot, aes_string(x = xvar, y = yvar, fill = fillvar), subtitle = subtitle) + 
        geom_raster() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        facet_grid(as.formula(paste(facetx, "~", facety)), labeller = label_both) +
        labs(caption = mycaption) +
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
        filename <- paste0(DateTimeStamp, "output", fillvarname, iteration, ".png")
        if(file.exists(filename)) {
          warning("File already exists, not saved again!")
        } else {
          ggsave(filename)
        }
      }
      iteration <- iteration + 1
    }
  }
}

RunOverTime <- function(parms = Mydf, verbose = FALSE, type = "Pair", ...) {
  state <- c(Nutr = parms[["NutrInit"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumInit"]],
             TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0, MttLum = 0,
             DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]],
             TransWall = 0, MdrWall = 0, MdtWall = 0, MrtWall = 0, MttWall = 0)
  out2 <- ode(t = times, y = state, func = ModelPairsNutr, parms = parms, verbose = verbose)
  out2 <- cbind(out2, TotalDLum = NA, TotalRLum = NA, TotalTransLum = NA,
                TotalDWall = NA, TotalRWall = NA, TotalTransWall = NA)
  out2[, "TotalDLum"] <- out2[, "DLum"] + out2[, "MdrLum"] + out2[, "MdtLum"]
  out2[, "TotalRLum"] <- out2[, "RLum"] + out2[, "MdrLum"] + out2[, "MrtLum"]
  out2[, "TotalTransLum"] <- out2[, "TransLum"] + out2[, "MdtLum"] + out2[, "MrtLum"] + 2*out2[, "MttLum"]
  out2[, "TotalDWall"] <- out2[, "DWall"] + out2[, "MdrWall"] + out2[, "MdtWall"]
  out2[, "TotalRWall"] <- out2[, "RWall"] + out2[, "MdrWall"] + out2[, "MrtWall"]
  out2[, "TotalTransWall"] <- out2[, "TransWall"] + out2[, "MdtWall"] + out2[, "MrtWall"] + 2*out2[, "MttWall"]
  
  EqAfterInvasion <- tail(out2, 1)
  if(verbose == TRUE) {
    print(diagnostics(out2))
    print(attributes(out2))
  }
  PlotOverTime(plotdata = out2, parms = parms, type = type, verbose = verbose, saveplot = saveplots)
  stateBulk <- c(Nutr = parms[["NutrInit"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumInit"]], TransLum = 0,
                 DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]], TransWall = 0)
  out2bulk <- ode(t = times, y = stateBulk, func = ModelBulkNutr, parms = parms, verbose = verbose)
  EqAfterInvasionBulk <- tail(out2bulk, 1)
  if(verbose == TRUE) {
    print(diagnostics(out2bulk))
    print(attributes(out2bulk))
  }
  PlotOverTime(plotdata = out2bulk, parms = parms, type = "Bulk", saveplot = saveplots)
  EqAfterInvasionTotal <- cbind(EqAfterInvasion, EqAfterInvasionBulk)
  names(EqAfterInvasionTotal) <- c(colnames(EqAfterInvasion),
                                   paste0(colnames(EqAfterInvasionBulk), "Bulk"))
  return(EqAfterInvasionTotal)
}

PlotOverTime <- function(plotdata = out2, parms = parms, type = "Pair", saveplot = saveplots, ...) {
  subtitle <- paste0("kp=", parms[["kp"]], ", ", parms[["kpWall"]],
                     " kn=", parms[["kn"]], ", ", parms[["knWall"]],
                     " gd=", parms[["gd"]], " gt=", parms[["gt"]],
                     " gdbulk=", signif(parms[["gdbulkLum"]], 3),
                     ", ", signif(parms[["gdbulkWall"]], 3), 
                     " gtbulk=", signif(parms[["gtbulkLum"]], 3),
                     ", ", signif(parms[["gtbulkWall"]], 3),
                     " cd=", parms[["cd"]], " ct=", parms[["ct"]],
                     " bR=", parms[["bR"]], " NI=", parms[["NI"]],
                     " NutrConv=", parms[["NutrConv"]], " w=", parms[["w"]]
  )
  if(type == "Pair") {
    maintitle <- "Pair model"
    mycol <- mycolpairs
    mylty <- myltypairs
    mylwd <- c(rep(2, 8), rep(1, 7))
    SelectedColumns <- c("Nutr", "DLum", "RLum", "TransLum", "MdrLum", "MdtLum", "MrtLum", "MttLum",
                         "DWall", "RWall", "TransWall", "MdrWall", "MdtWall", "MrtWall", "MttWall")
  } else {
    mycol <- mycolother
    mylty <- myltyother
    mylwd <- c(rep(2, 4), rep(1, 3))
    if(type == "Total") {
      maintitle <- "Pair model, totals"
      SelectedColumns <- c("Nutr", "TotalDLum", "TotalRLum", "TotalTransLum",
                           "TotalDWall", "TotalRWall", "TotalTransWall")
    } else {
      maintitle <- "Bulk model"
      SelectedColumns <- NULL
    }
  }
  if(saveplot == TRUE) {
    filename <- paste0(DateTimeStamp, "output", gsub(" ", "", maintitle), ".png")
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      png(filename = filename)
    }
    matplot.deSolve(plotdata, which = SelectedColumns, main = maintitle,
                    sub = subtitle, ylim = myylim, log = if(yaxislog == 1) {"y"},
                    col = mycol, lty = mylty, lwd = mylwd,
                    legend = list(x = "bottomright"))
    grid()
    if(file.exists(filename) == FALSE) {
      dev.off()
    }
  } else {
    matplot.deSolve(plotdata, which = SelectedColumns, main = maintitle,
                    sub = subtitle, ylim = myylim, log = if(yaxislog == 1) {"y"},
                    col = mycol, lty = mylty, lwd = mylwd,
                    legend = list(x = "bottomright"))
    grid()
  }
}

#### Main script ####
CheckParms <- c(NILum = NILumSet, NIWall = NIWallSet, wLum = wLumSet, wWall = wWallSet, 
                NutrConv = NutrConvSet, bR = bRSet,
                MigrLumWall = MigrLumWallSet, MigrWallLum = MigrWallLumSet, ScaleAreaPerVol = ScaleAreaPerVolSet,
                DInitLum = DInitLumSet, DInitWall = DInitWallSet, 
                cd = cdSet, ct = ctSet,
                kp = kpSet, kpWall = kpWallSet, kn = knSet, knWall = knWallSet,
                gd = gdSet, gt = gtSet)
warntext <- paste("Parameter(s)", names(which(CheckParms <= 0)), "contain(s) non-positive values.")
if(any(CheckParms <= 0)) warning(warntext)
if(any(c(cdSet, ctSet) >= 1)) warning("Costs should be larger than 0 and smaller than 1.")

TotalIterations <- length(NILumSet)*length(NIWallSet)*length(wLumSet)*length(wWallSet)*
  length(NutrConvSet)*length(bRSet)*
  length(MigrLumWallSet)*length(MigrWallLumSet)*length(ScaleAreaPerVolSet)*
  length(DInitLumSet)*length(DInitWallSet)*
  length(cdSet)*length(ctSet)*
  length(kpSet)*length(kpWallSet)*length(knSet)*length(knWallSet)*
  length(gdSet)*length(gtSet)
TotalIterations

## Determine plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(NILum = NILumSet, NIWall = NIWallSet, wLum = wLumSet, wWall = wWallSet,
                      NutrConv = NutrConvSet, bR = bRSet, 
                      MigrLumWall = MigrLumWallSet, MigrWallLum = MigrWallLumSet,
                      ScaleAreaPerVol = ScaleAreaPerVolSet)
dim(MyData)

# Add equilibrium values of the one-compartment model as initial state values to
# run to the plasmid-free equilibrium
MyData <-  mutate(MyData, 
                  NutrLumGuess = wLum/bR,
                  RLumGuess = (NILum - wLum/bR)/NutrConv,
                  NutrWallGuess = wWall/bR,
                  RWallGuess = (NIWall - wWall/bR)/NutrConv
)

if(any(select(MyData, ends_with("Guess")) < 0)) {
  warning("Some of the initial state values are negative!")
}

# Determine the cell-containing plasmid-free equilibrium (NutrLum*, RLum*,
# NutrWall*, RWall*) by running the model with only nutrients and recipients to
# the plasmid-free equilibrium.
Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = RunToPlamidfreeEq))
if(any(Eqplasmidfree[, "steadyInit"] != 1)) {
  warning("Steady-state has not been reached for all plasmid-free equilibria!")
}

if(any(Eqplasmidfree <= 0)) {
  RowsNegativeEq <- unique(which(Eqplasmidfree <= 0, arr.ind = TRUE)[, 1])
  ColnamesNegativeEq <- colnames(Eqplasmidfree)[unique(which(Eqplasmidfree <= 0, arr.ind = TRUE)[, 2])]
  warning("Plasmid-free equilibrium contains non-positive values in columns '", paste(ColnamesNegativeEq, collapse = "' and '"),
          "'.\nThe data will be included in the calculations anyway!
  This concerns the following rows of the dataframe: ", paste(RowsNegativeEq, collapse = ", "))
}

MyData <- cbind(MyData, Eqplasmidfree)
print("Plasmid-free equilibrium determined:")
print(Sys.time())

### Checking if the biomass in the two compartments are comparable 
summary(MyData[, "RLumInit"] / MyData[, "RWallInit"])
summary(MyData[, "NutrLumInit"] / MyData[, "NutrWallInit"])

limitsREq <- log10(c(min(c(MyData$RLumInit, MyData$RWallInit)), max(c(MyData$RLumInit, MyData$RWallInit))))
limitsNutrEq <- log10(c(min(c(MyData$NutrLumInit, MyData$NutrWallInit)), max(c(MyData$NutrLumInit, MyData$NutrWallInit))))

DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")
CreatePlot(fillvar = "log10(RLumInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsREq)
CreatePlot(fillvar = "log10(RWallInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsREq)
CreatePlot(fillvar = "log10(NutrLumInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsNutrEq)
CreatePlot(fillvar = "log10(NutrWallInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsNutrEq)

CreatePlot(fillvar = "RLumInit/RWallInit", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall")
CreatePlot(fillvar = "NutrLumInit/NutrWallInit", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall")

## Approximate gdbulk and gtbulk in the lumen
MyData <- expand_grid(MyData, DInitLum = DInitLumSet, DInitWall = DInitWallSet,
                      kp = kpSet, kn = knSet, gd = gdSet, gt = gtSet)
DataEstConjBulk <- t(apply(X = MyData, MARGIN = 1, FUN = EstConjBulkLum))

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulkLum <- unname(MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor))

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] + DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] + DataEstConjBulk[, "TransMrt"]
gtbulkLum <- unname(MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans))
MyData <- cbind(MyData, gdbulkLum = gdbulkLum, gtbulkLum = gtbulkLum)

## Approximate gdbulk and gtbulk at the wall
MyData <- expand_grid(MyData, kpWall = kpWallSet, knWall = knWallSet)

# Replace columns kn and kn with the values in the columns kpWall and knWall in
# new dataframe to be used to estimate bulk conjugation rates at the wall
MyDataWall <- select(MyData, !c(kp, kn, gdbulkLum, gtbulkLum))
MyDataWall <- rename(MyDataWall, kp = kpWall, kn = knWall)
DataEstConjBulk <- t(apply(X = MyDataWall, MARGIN = 1, FUN = EstConjBulkWall))
DataEstConjBulk <- as.data.frame(DataEstConjBulk)

DataEstConjBulk <- mutate(DataEstConjBulk,
                                 TotalDEstConjBulkDonor = DonorD + DonorMdr + DonorMdt,
                                 TotalREstConjBulkDonor = DonorR + DonorMdr + DonorMrt,
                                 gdbulkWallpart = DonorMdr / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor))
gdbulkWall <- unname(MyData[, "gd"] * DataEstConjBulk[, "gdbulkWallpart"])

DataEstConjBulk <- mutate(DataEstConjBulk,
                          TotalTransEstConjBulkTrans = TransTrans + TransMrt + 2*TransMtt,
                          TotalREstConjBulkTrans = TransR + TransMrt,
                          gtbulkWallpart = TransMrt / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans))
gtbulkWall <- unname(MyData[, "gt"] * DataEstConjBulk[, "gtbulkWallpart"])

MyData <- cbind(MyData, gdbulkWall = gdbulkWall, gtbulkWall = gtbulkWall)

print("Bulk-conjugation rates estimated:")
print(Sys.time())

#### The remainder of the script still has to be adjusted to the new model ####


## WERKT NIET ?
CumRowIndex <- NULL
iteration <- 1
dataplottotal <- MyData
for(MigrLumWallSubset in unique(dataplottotal[, "MigrLumWall"])) {
  subtitle <- paste0("MigrLumWall= ", MigrLumWallSubset)
  RowIndex <- dataplottotal[, "MigrLumWall"] == MigrLumWallSubset
  dataplot <- dataplottotal[RowIndex, ]
  CreatePlot2(fillvar = "SignDomEigVal", gradient2 = 1, limits = c(-1, 1), facetx = "kpWall", facety = "knWall",
              save = TRUE)
  iteration <- iteration + 1
}


# calculate (or approximate?) eigenvalues
MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)

MyInfoEigVal <- t(apply(MyData, MARGIN = 1, FUN = CalcEigenvalues))
MyData <- cbind(MyData, MyInfoEigVal)

print("Eigenvalues estimated:")
print(Sys.time())

write.csv(MyData, file = paste0(DateTimeStamp, "outputnosimulationtwocompartment.csv"),
          quote = FALSE, row.names = FALSE)

CreatePlot(fillvar = "SignDomEigVal", gradient2 = 1, limits = c(-1, 1), facetx = "MigrLumWall", facety = ".", save = TRUE)

CreatePlot2(fillvar = "SignDomEigVal", gradient2 = 1, xvar = "MigrLumWall", yvar = "w")
CreatePlot2(fillvar = "SignDomEigVal", gradient2 = 1, xvar = "MigrLumWall", yvar = "MigrWallLum")

limitsREq <- log10(c(min(c(MyData$RLumInit, MyData$RWallInit)), max(c(MyData$RLumInit, MyData$RWallInit))))

CreatePlot2(fillvar = "log10(RLumInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", limits = limitsREq, save = TRUE)
CreatePlot2(fillvar = "log10(RLumInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", limits = limitsREq, save = TRUE)
CreatePlot2(fillvar = "NutrInit", xvar = "MigrLumWall", yvar = "MigrWallLum", save = TRUE)

### Some plotting ###
# See the one-compartment script for more plots
CreatePlot(fillvar = "SignDomEigVal", gradient2 = 1, limits = c(-1, 1))
CreatePlot(fillvar = "SignDomEigValBulk", gradient2 = 1, limits = c(-1, 1))

DomEigVals <- c(MyData$DomEigVal, MyData$DomEigValBulk)
limitseigenvalues <- log10(range(DomEigVals[DomEigVals > 0]))
limitseigenvalues <- c(floor(limitseigenvalues[1]), ceiling(limitseigenvalues[2]))

CreatePlot(fillvar = "log10(DomEigVal)", limits = limitseigenvalues)
CreatePlot(fillvar = "log10(DomEigValBulk)", limits = limitseigenvalues)
CreatePlot(fillvar = "DomEigVal/DomEigValBulk")
CreatePlot(fillvar = "DomEigVal-DomEigValBulk")

limitsbulkrates <- c(floor(min(log10(c(MyData$gdbulkLum, MyData$gtbulkLum, MyData$gdbulkWall, MyData$gtbulkWall)))),
                     ceiling(max(log10(c(MyData$gdbulkLum, MyData$gtbulkLum, MyData$gdbulkWall, MyData$gtbulkWall)))))
CreatePlot(fillvar = "log10(gdbulkLum)", limits = limitsbulkrates)
CreatePlot(fillvar = "log10(gtbulkLum)", limits = limitsbulkrates)
CreatePlot(fillvar = "log10(gdbulkWall)", limits = limitsbulkrates)
CreatePlot(fillvar = "log10(gtbulkWall)", limits = limitsbulkrates)

ggplot(data = MyData, aes(x = log10(kp), y = log10(kn), fill = factor(SignDomEigVal))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(cd + gd ~ ct + gt, labeller = label_both) +
  labs(caption = DateTimeStamp) +
  theme(legend.position = "bottom", plot.caption = element_text(vjust = 20)) +
  scale_fill_manual(values = c("-1" = "darkblue", "1" = "darkred"),
                    name = "Plasmid-free equilibrium",
                    labels = c("stable", "unstable"))

# If invasion is possible, run simulation to see how many bacteria of each
# population are present at equilibrium

IndexSimulation <- which(MyData$SignDomEigVal != -1)
print(paste(length(IndexSimulation), "simulations to run for the pair-formation model"))
ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))
InputSimulationPairs <- MyData[IndexSimulation, ColumnsToSelect]
OutputSimulationPairs <- t(apply(X = InputSimulationPairs, MARGIN = 1,
                                 FUN = SimulationPairs))

if(length(IndexSimulation) < nrow(MyData)) {
  NoSimulationNeeded <- cbind(time = 0, steady = 1, Nutr = MyData[-IndexSimulation, "NutrInit"],
                              DLum = 0, RLum = MyData[-IndexSimulation, "RLumInit"],
                              TransLum = 0,  MdrLum = 0, MdtLum = 0, MrtLum = 0,
                              MttLum = 0, DWall = 0,
                              RWall = MyData[-IndexSimulation, "RLumInit"],
                              TransWall = 0, MdrWall = 0, MdtWall = 0,
                              MrtWall = 0, MttWall = 0)
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
  NoSimulationNeededBulk <- cbind(timeBulk = 0, steadyBulk = 1,
                                  NutrBulk = MyData[-IndexSimulation, "NutrInit"],
                                  DLumBulk = 0,
                                  RLumBulk = MyData[-IndexSimulation, "RLumInit"],
                                  TransLumBulk = 0, DWallBulk = 0,
                                  RWallBulk = MyData[-IndexSimulation, "RLumInit"],
                                  TransWallBulk = 0)
  MyData <- rbind(cbind(MyData[IndexSimulationBulk, ], OutputSimulationBulk),
                  cbind(MyData[-IndexSimulationBulk, ], NoSimulationNeededBulk))
} else {
  MyData <- cbind(MyData, OutputSimulationBulk)
}
if(any(MyData$steadyBulk == 0)) warning("Steady-state has not always been reached")

print("Bulk-conjugation model completed running:")
print(Sys.time())

MyData <- cbind(MyData, TotalDLum = NA, TotalRLum = NA, TotalTransLum = NA, TotalPlasmidLum = NA, TotalBioLum = NA,
                TotalDWall = NA, TotalRWall = NA, TotalTransWall = NA, TotalPlasmidWall = NA, TotalBioWall = NA,
                TotalPlasmidLumBulk = NA, TotalBioLumBulk = NA, TotalPlasmidWallBulk = NA, TotalBioWallBulk = NA)
MyData[, "TotalDLum"] <- MyData[, "DLum"] + MyData[, "MdrLum"] + MyData[, "MdtLum"]
MyData[, "TotalRLum"] <- MyData[, "RLum"] + MyData[, "MdrLum"] + MyData[, "MrtLum"]
MyData[, "TotalTransLum"] <- MyData[, "TransLum"] + MyData[, "MdtLum"] + MyData[, "MrtLum"] + 2*MyData[, "MttLum"]
MyData[, "TotalPlasmidLum"] <- MyData[, "TotalDLum"] + MyData[, "TotalTransLum"]
MyData[, "TotalBioLum"] <- MyData[, "TotalRLum"] + MyData[, "TotalPlasmidLum"]
MyData[, "TotalDWall"] <- MyData[, "DWall"] + MyData[, "MdrWall"] + MyData[, "MdtWall"]
MyData[, "TotalRWall"] <- MyData[, "RWall"] + MyData[, "MdrWall"] + MyData[, "MrtWall"]
MyData[, "TotalTransWall"] <- MyData[, "TransWall"] + MyData[, "MdtWall"] + MyData[, "MrtWall"] + 2*MyData[, "MttWall"]
MyData[, "TotalPlasmidWall"] <- MyData[, "TotalDWall"] + MyData[, "TotalTransWall"]
MyData[, "TotalBioWall"] <- MyData[, "TotalRWall"] + MyData[, "TotalPlasmidWall"]
MyData[, "TotalPlasmidLumBulk"] <- MyData[, "DLumBulk"] + MyData[, "TransLumBulk"]
MyData[, "TotalBioLumBulk"] <- MyData[, "RLumBulk"] + MyData[, "TotalPlasmidLumBulk"]
MyData[, "TotalPlasmidWallBulk"] <- MyData[, "DWallBulk"] + MyData[, "TransWallBulk"]
MyData[, "TotalBioWallBulk"] <- MyData[, "RWallBulk"] + MyData[, "TotalPlasmidWallBulk"]

write.csv(MyData, file = paste0(DateTimeStamp, "outputtwocompartment.csv"),
          quote = FALSE, row.names = FALSE)

CreatePlot(fillvar = "TotalRLum/TotalBioLum")
CreatePlot(fillvar = "TotalRWall/TotalBioWall")

CreatePlot(fillvar = "RLumBulk/TotalBioLumBulk")
CreatePlot(fillvar = "RWallBulk/TotalBioWallBulk")

CreatePlot(fillvar = "(RBulk/TotalBioBulk) / (TotalR/TotalBio)")

##### Create plots over time #####
myltypairs <- c(lty = c(3, rep(c(1, 2, 1, 1, 1, 1, 1), 2)))
myltyother <- c(lty = c(3, rep(c(1, 2, 1), 2)))
mycolpairs <- c("black", rep(c("purple", "green1", "red", "yellow", "hotpink", "blue", "cyan"), 2))
mycolother <- c("black", rep(c("purple", "green1", "red"), 2))
myylim <- c(1E-6, 1E7)
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
verbose <- 0 # if verbose == 1, diagnositics on the simulations are printed
Mytmax <- c(500)
Mytstep <- c(0.1)
TheseRows <- 1:nrow(MyData) # Rows to use for simulations over time

Mydf <- MyData[TheseRows, ColumnsToSelect]
TotalIterations <- length(TheseRows)
print(TotalIterations)

# Times for which output of the simulation is wanted. Note that the used
# ode-solvers are variable-step methods, so the times in times are NOT the only
# times at which integration is performed. See help(diagnostics.deSolve()) and
# help(lsodar()) for details.
times <- seq(from = 0, to = Mytmax, by = Mytstep)

# To see the dynamics of the different populations
EqAfterInvasion <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime, type = "Pair"))

# To compare total numbers of donors, recipients, and transconjugants in the
# output of the pair-formation model with the bulk-formation model
EqAfterInvasion <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime, type = "Total"))
