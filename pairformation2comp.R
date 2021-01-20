#### Two-compartment pair-formation model of conjugation ####

#### To do ####
# Also see the 'To do' section in the pair-formation script for the one-compartment model

# Instead of checking if the plasmid-free equilibrium that has been found is valid,
# I could derive the invasion-criterion of the cell-free equilibrium
# (i.e. Nutr* = NI, D*, R*, ... = 0) to determine if cells can exist.

# ScaleAreaPerVol is now used to check input and written to the .csv-file, that
# has to be replaced by VLum and VWall.

# Can I use the migration rates to set better states for the plasmid-free equilibrium

# See Macken 1994 'The dynamics of bacteria-plasmid systems' for analytic treatment
# and a way to get rid of the nutrient equation (but is that feasible if I want
# to modify nutrient levels over time?)

# In filtering data, using filter(MyData, kpWall == i, knWall == j) is error-prone
# because comparing two floating-point vectors gives troubles.
# Use dplyr::near(kpWall, i), or maybe base-R isTRUE(all.equal(...)).

# Check if using MigrLumWall = MigrWallLum = 0 en DInitWall = 0 leads to same 
# results as the single-compartment model.

# Also use ModelEstConjBulkDonor to approximate bulk conjugation rates from the transconjugant
# (note: than Mdr should again be added in calculation of TotalREstConjBulkTrans)

# The calculations of TotalDEstConjBulkDonor, TotalREstConjBulkDonor, TotalTransEstConjBulkDonor
# ect. for approximating bulk-rates can be moved inside the EstConjBulk function ?

# In calculation of bulkrates use knLum, knWall instead of kn and knWall, then
# use kn = knLum as function argument for lumen and kn = knWall as function argument for wall.
# This will prevent creating 'DataWall' as separate dataframe

# Saving plots in PlotOverTime() (called from RunOverTime()) fails after the
# first plot because names are the same. Add 'Iteration' as column to data
# and append the value from that column to the name? That would als make it
# easier to see which plot is which iteration (so also print it in graph subtitles).

# Now a loop is used to get a dataframe with percentage and counts of parameter
# combinations for which the plasmid can, or cannot, invade. That could go into
# a function. 

# dev.off() is not always called when needed
 
# Total plasmid and biomass are not calculated when simulating over time

# An alternative to using the labelling function and supplying names would be to
# rename the variables before plotting. This has as drawback that names cannot
# contain spaces, and I don't know how automated wrapping long names behaves.
# MyData <- rename(MyData, "Nutr_conc_inflow" = "NI", "Washout_rate" = "w")

## To check if other parameters than kp, kn, gd, gt, cd, ct which are separately
# shown in the plot have multiple values (which leads to plotting multiple values
# on top of each other) the following can be used:
# # Parms <- list(DInitLumSet = DInitLumSet, DInitWallSet = DInitWallSet, bRSet = bRSet,
# NISet = NISet, NutrConvSet = NutrConvSet, MigrLumWallSet = MigrLumWallSet,
# MigrWallLumSet = MigrWallLumSet,
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
library(deSolve) # Integrate differential equations with results over time.
library(ggplot2) # To create plots
library(RColorBrewer) # To get better color schemes
library(rootSolve) # Integration, obtaining jacobian matrix and eigenvalues.
library(tidyr) # expand.grid() with dataframe as input
library(dplyr) # mutate(), filter(), near()

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
MyColorBrew <- rev(brewer.pal(11, "Spectral")) # examples: display.brewer.all()

#### Functions ####
# Model containing only nutrients and recipients.
ModelBulkNutrPlasmidfree <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutrLum <- ((NILum - NutrLum)*wLum - NutrConv*bR*NutrLum*RLum/(Ks + NutrLum))*VLum
    dRLum <- (bR*NutrLum/(Ks + NutrLum) - wLum - MigrLumWall)*RLum*VLum + MigrWallLum*RWall*VWall
    dNutrWall <- ((NIWall - NutrWall)*wNutrWall - NutrConv*bR*NutrWall*RWall/(Ks + NutrWall))*VWall
    dRWall <- (bR*NutrWall/(Ks + NutrWall) - wWall - MigrWallLum)*RWall*VWall + MigrLumWall*RLum*VLum
    return(list(c(dNutrLum, dRLum, dNutrWall, dRWall)))
  })
}

# Model to determine plasmid-free equilibrium, expressed as milligram nutrients
# and number of recipient bacteria in the lumen and wall compartments.
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
# are not included in this model. The state variables are the number of cells in
# the lumen or wall compartment, not the concentration (i.e., units are cells,
# not cells per mL).
ModelEstConjBulkDonor <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dD <- (- kp*D*R + kn*(Mdr + Mdt))*VComp
    dR <- (- kp*R*(D + Trans) + kn*(Mdr + Mrt))*VComp
    dTrans <- (- kp*R*Trans + kn*(Mdt + Mrt))*VComp
    dMdr <- (kp*D*R - kn*Mdr - gd*Mdr)*VComp
    dMdt <- (gd*Mdr - kn*Mdt)*VComp
    dMrt <- (kp*R*Trans - kn*Mrt)*VComp
    return(list(c(dD, dR, dTrans, dMdr, dMdt, dMrt)))
  })
}

# ODE-model used to approximate the bulk-conjugation rate of the transconjugant.
# Nutrients, growth, washout, and donors are not included in this model. The
# state variables are the number of cells in the lumen or wall compartment, not
# the concentration (i.e., units are cells, not cells per mL).
ModelEstConjBulkTrans <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR <- (- kp*R*Trans + kn*Mrt)*VComp
    dTrans <- (- kp*R*Trans + kn*(Mrt + 2*Mtt))*VComp
    dMrt <- (kp*R*Trans - kn*Mrt - gt*Mrt)*VComp
    dMtt <- (gt*Mrt - kn*Mtt)*VComp
    return(list(c(dR, dTrans, dMrt, dMtt)))
  })
}

# Functions to estimate bulk-conjugation rates by running simulations with the
# adjusted pair-formation models for a short time (i.e., tail(timesEstConj, 1)
# hours) and calculate approximations of gdbulk and gtbulk from the output,
# following Zhong's approach for the calculations. The number of donors is
# supplied as number of cells per mL so should be multiplied by the appropriate
# volumes to be used as state. RLumInit and RWallInit have been calculated as
# the number of cells, so do not need to be multiplied by volume again.
# 
# NOTE: now I actually calculate the bulk-rates as I used to do (so using densities
# and concentrations), and only convert
# them when using them in the bulk-conjugation models.
# 
EstConjBulkLum <- function(MyData) {
  with(as.list(MyData), {
    if(MyData[["DInitLum"]] == 0) {
      state <- c(D = MyData[["DInitWall"]], R = MyData[["RLumInit"]]/VLum,
                 Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    } else {
      state <- c(D = MyData[["DInitLum"]], R = MyData[["RLumInit"]]/VLum,
                 Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    }
    
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    if(MyData[["DInitLum"]] == 0) {
      state <- c(R = MyData[["RLumInit"]]/VLum, Trans = MyData[["DInitWall"]],
                 Mrt = 0, Mtt = 0)
    } else {
      state <- c(R = MyData[["RLumInit"]]/VLum, Trans = MyData[["DInitLum"]],
                 Mrt = 0, Mtt = 0)
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
      state <- c(D = MyData[["DInitLum"]], R = MyData[["RWallInit"]]/VWall,
                 Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    } else {
      state <- c(D = MyData[["DInitWall"]], R = MyData[["RWallInit"]]/VWall,
                 Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    }
    
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    if(MyData[["DInitWall"]] == 0) {
      state <- c(R = MyData[["RWallInit"]]/VWall, Trans = MyData[["DInitLum"]],
                 Mrt = 0, Mtt = 0)
    } else {
      state <- c(R = MyData[["RWallInit"]]/VWall, Trans = MyData[["DInitWall"]],
                 Mrt = 0, Mtt = 0)
    }
    
    DataEstConjBulkTrans <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)
    
    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
  })
}

## NOTE: have to  CHECK IF SCALING for migration in two-comp model DID NOT INTRODUCE ERRORS


# Bulk-conjugation model, with inflow, outflow, conversion by bacteria for
# nutrient equations, and growth, washout, migration from and to the lumen and
# wall compartment, conjugation from donors and transconjugants for bacterial
# equations. State variables are expressed as milligram nutrients and number of
# bacteria.
ModelBulkNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutrLum <- ((NILum - NutrLum)*wLum - NutrConv*bR*NutrLum*((1 - cd)*DLum + RLum + (1 - ct)*TransLum)/(Ks + NutrLum))*VLum
    dDLum <- ((1 - cd)*bR*NutrLum/(Ks + NutrLum) - wLum - MigrLumWall)*DLum*VLum + MigrWallLum*DWall*VWall
    dRLum <- (bR*NutrLum/(Ks + NutrLum) - wLum - MigrLumWall)*RLum*VLum + MigrWallLum*RWall*VWall -
      (gdbulkLum*DLum + gtbulkLum*TransLum)*RLum*VLum
    dTransLum <- ((1 - ct)*bR*NutrLum/(Ks + NutrLum) - wLum - MigrLumWall)*TransLum*VLum + MigrWallLum*TransWall*VWall +
      (gdbulkLum*DLum + gtbulkLum*TransLum)*RLum*VLum

    dNutrWall <- ((NIWall - NutrWall)*wNutrWall - NutrConv*bR*NutrWall*((1 - cd)*DWall + RWall + (1 - ct)*TransWall)/(Ks + NutrWall))*VWall
    dDWall <- ((1 - cd)*bR*NutrWall/(Ks + NutrWall) - wWall - MigrWallLum)*DWall*VWall + MigrLumWall*DLum*VLum
    dRWall <- (bR*NutrWall/(Ks + NutrWall) - wWall - MigrWallLum)*RWall*VWall + MigrLumWall*RLum*VLum -
      (gdbulkWall*DWall + gtbulkWall*TransWall)*RWall*VWall
    dTransWall <- ((1 - ct)*bR*NutrWall/(Ks + NutrWall) - wWall - MigrWallLum)*TransWall*VWall + MigrLumWall*TransLum*VLum +
      (gdbulkWall*DWall + gtbulkWall*TransWall)*RWall*VWall
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
# (Zhong 2010). The structure of the two compartments and migration between them
# is based on Imran (2005), but I added nutrient inflow and washout of nutrients
# and bacteria from the wall compartment. I also added costs in growth for
# plasmid-bearing bacteria. State variables are expressed as milligram nutrients and number of
# bacteria.

## NOTE: have to  CHECK IF SCALING for migration in two-comp model DID NOT INTRODUCE ERRORS

ModelPairsNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutrLum <- ((NILum - NutrLum)*wLum -
      NutrConv*bR*NutrLum*(
        (1 - cd)*(DLum + MdrLum + MdtLum) + (RLum + MdrLum + MrtLum) + (1 - ct)*(TransLum + MdtLum + MrtLum + 2*MttLum)/(Ks + NutrLum)
      ))*VLum
    dDLum <- ((1 - cd)*bR*NutrLum*(DLum + MdrLum + MdtLum)/(Ks + NutrLum) - (wLum + MigrLumWall + kp*RLum)*DLum +
      kn*(MdrLum + MdtLum))*VLum + MigrWallLum*DWall*VWall
    dRLum <- (bR*NutrLum*(RLum + MdrLum + MrtLum)/(Ks + NutrLum) - (wLum + MigrLumWall + kp*(DLum + TransLum))*RLum +
      kn*(MdrLum + MrtLum))*VLum + MigrWallLum*RWall*VWall
    dTransLum <- ((1 - ct)*bR*NutrLum*(TransLum + MdtLum + MrtLum + 2*MttLum)/(Ks + NutrLum) - (wLum + MigrLumWall + kp*RLum)*TransLum +
      kn*(MdtLum + MrtLum + 2*MttLum))*VLum + MigrWallLum*TransWall*VWall
    dMdrLum <- (kp*DLum*RLum - (kn + gd + wLum + MigrLumWall)*MdrLum)*VLum + MigrWallLum*MdrWall*VWall
    dMdtLum <- (gd*MdrLum - (kn + wLum + MigrLumWall)*MdtLum)*VLum + MigrWallLum*MdtWall*VWall
    dMrtLum <- (kp*RLum*TransLum - (kn + gt + wLum + MigrLumWall)*MrtLum)*VLum + MigrWallLum*MrtWall*VWall
    dMttLum <- (gt*MrtLum - (kn + wLum + MigrLumWall)*MttLum)*VLum + MigrWallLum*MttWall*VWall
    
    dNutrWall <- ((NIWall - NutrWall)*wNutrWall -
      NutrConv*bR*NutrWall*(
        (1 - cd)*(DWall + MdrWall + MdtWall) + (RWall + MdrWall + MrtWall) + (1 - ct)*(TransWall + MdtWall + MrtWall + 2*MttWall)/(Ks + NutrWall)
      ))*VWall
    dDWall <- ((1 - cd)*bR*NutrWall*(DWall + MdrWall + MdtWall)/(Ks + NutrWall) - (wWall + MigrWallLum + kpWall*RWall)*DWall +
      knWall*(MdrWall + MdtWall))*VWall + MigrLumWall*DLum*VLum
    dRWall <- (bR*NutrWall*(RWall + MdrWall + MrtWall)/(Ks + NutrWall) - (wWall + MigrWallLum + kpWall*(DWall + TransWall))*RWall +
      knWall*(MdrWall + MrtWall))*VWall + MigrLumWall*RLum*VLum
    dTransWall <- ((1 - ct)*bR*NutrWall*(TransWall + MdtWall + MrtWall + 2*MttWall)/(Ks + NutrWall) - (wWall + MigrWallLum + kpWall*RWall)*TransWall +
      knWall*(MdtWall + MrtWall + 2*MttWall))*VWall + MigrLumWall*TransLum*VLum
    dMdrWall <- (kpWall*DWall*RWall - (knWall + gd + wWall + MigrWallLum)*MdrWall)*VWall + MigrLumWall*MdrLum*VLum
    dMdtWall <- (gd*MdrWall - (knWall + wWall + MigrWallLum)*MdtWall)*VWall + MigrLumWall*MdtLum*VLum
    dMrtWall <- (kpWall*RWall*TransWall - (knWall + gt + wWall + MigrWallLum)*MrtWall)*VWall + MigrLumWall*MrtLum*VLum
    dMttWall <- (gt*MrtWall - (knWall + wWall + MigrWallLum)*MttWall)*VWall + MigrLumWall*MttLum*VLum
    
    return(list(c(dNutrLum, dDLum, dRLum, dTransLum, dMdrLum, dMdtLum, dMrtLum, dMttLum,
                  dNutrWall, dDWall, dRWall, dTransWall, dMdrWall, dMdtWall, dMrtWall, dMttWall)))
  })
}

# Numerically estimate the Jacobian matrix of the plasmid-free equilibrium of
# the models, then approximate the eigenvalues of this matrix.
# The maximum real part of the eigenvalues is used to determine stability.
CalcEigenvalues <- function(MyData) {
  parms <- MyData
  EqFull <- c(NutrLum = MyData[["NutrLumInit"]], DLum = 0, RLum = MyData[["RLumInit"]],
              TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0, MttLum = 0,
              NutrWall = MyData[["NutrWallInit"]], DWall = 0, RWall = MyData[["RWallInit"]],
              TransWall = 0, MdrWall = 0, MdtWall = 0, MrtWall = 0, MttWall = 0)
  EigVal <- eigen(x = jacobian.full(y = EqFull, func = ModelPairsNutr, parms = parms),
                  symmetric = FALSE, only.values = TRUE)$values
  ComplexEigVal <- is.complex(EigVal) 
  EigVal <- Re(EigVal)
  names(EigVal) <- paste0("Eigval", 1:length(EigVal))
  DomEigVal <- max(EigVal)
  SignDomEigVal <- sign(DomEigVal)
  SignEigValEqual <- identical(rep(SignDomEigVal, length(EigVal)), sign(Re(EigVal)))
  
  EqFullBulk <- c(NutrLum = MyData[["NutrLumInit"]], DLum = 0, RLum = MyData[["RLumInit"]], TransLum = 0, 
                  NutrWall = MyData[["NutrWallInit"]], DWall = 0, RWall = MyData[["RWallInit"]], TransWall = 0)
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

# !! NOTE: Equations (and states?) HAVE NOT YET BEEN converted to milligram
# nutrients and number of cells by multiplication with the appropriate volumes !!

SimulationPairs <- function(InputSimulationPairs) {
  parms <- InputSimulationPairs
  state <- c(NutrLum = parms[["NutrLumInit"]], DLum = parms[["DInitLum"]],
             RLum = parms[["RLumInit"]], TransLum = 0, MdrLum = 0, MdtLum = 0,
             MrtLum = 0, MttLum = 0, NutrWall = parms[["NutrWallInit"]],
             DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]],
             TransWall = 0, MdrWall = 0, MdtWall = 0,
             MrtWall = 0, MttWall = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelPairsNutr,
                   parms = parms, stol = 1.25e-6, atol = atol)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

# !! NOTE: Equations (and states?) HAVE NOT YET BEEN converted to milligram
# nutrients and number of cells by multiplication with the appropriate volumes !!

SimulationBulk <- function(InputSimulationBulk, state) {
  parms <- InputSimulationBulk
  state <- c(NutrLum = parms[["NutrLumInit"]], DLum = parms[["DInitLum"]],
             RLum = parms[["RLumInit"]], TransLum = 0, NutrWall = parms[["NutrWallInit"]], 
             DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]], TransWall = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelBulkNutr,
                   parms = parms, stol = 2.5e-6, atol = atol)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

# Function to create and save heatmaps
CreatePlot <- function(dataplot = MyData, xvar = "log10(kp)", yvar = "log10(kn)",
                       fillvar = "factor(SignDomEigVal)",
                       filltype = "discrete", limits = NULL, 
                       labx = "Log10(attachment rate in lumen)",
                       laby = "Log10(detachment rate in lumen)",
                       filltitle, filllabels = c("No", "Yes"),
                       mytag = NULL,
                       manualvalues = c("TRUE" = "darkgreen", "FALSE" = "red"),
                       facetx = "gt + ct", facety = "gd + cd", as.table = TRUE,
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
      filename <- paste0(DateTimeStamp, "output", fillvarname, "twocomp", ".png")      
    }
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
        DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
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

# !! NOTE: Equations (and states?) HAVE NOT YET BEEN converted to milligram
# nutrients and number of cells by multiplication with the appropriate volumes !!

RunOverTime <- function(parms = Mydf, verbose = FALSE, type = "Pair", ...) {
  state <- c(NutrLum = parms[["NutrLumInit"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumInit"]],
             TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0, MttLum = 0,
             NutrWall = parms[["NutrWallInit"]],
             DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]],
             TransWall = 0, MdrWall = 0, MdtWall = 0, MrtWall = 0, MttWall = 0
             )
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
  stateBulk <- c(NutrLum = parms[["NutrLumInit"]], DLum = parms[["DInitLum"]],
                 RLum = parms[["RLumInit"]], TransLum = 0,
                 NutrWall = parms[["NutrWallInit"]], 
                 DWall = parms[["DInitWall"]], RWall = parms[["RWallInit"]], TransWall = 0)
  out2bulk <- ode(t = times, y = stateBulk, func = ModelBulkNutr, parms = parms, verbose = verbose)
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

PlotOverTime <- function(plotdata = out2, parms = parms, type = "Pair", saveplot = saveplots, ...) {
  subtitle <- paste0("kp=", parms[["kp"]], ", ", parms[["kpWall"]],
                     " kn=", parms[["kn"]], ", ", parms[["knWall"]],
                     " gd=", parms[["gd"]], " gt=", parms[["gt"]],
                     " gdbulk=", signif(parms[["gdbulkLum"]], 3),
                     ", ", signif(parms[["gdbulkWall"]], 3), 
                     " gtbulk=", signif(parms[["gtbulkLum"]], 3),
                     ", ", signif(parms[["gtbulkWall"]], 3),
                     " cd=", parms[["cd"]], " ct=", parms[["ct"]],
                     " bR=", parms[["bR"]], " NILum=", parms[["NILum"]],
                     " NIWall=", parms[["NIWall"]],
                     " NutrConv=", parms[["NutrConv"]],
                     " wLum=", parms[["wLum"]], " wWall=", parms[["wWall"]]
  )
  if(type == "Pair") {
    maintitle <- "Pair model"
    mycol <- mycolpairs
    mylty <- myltypairs
    mylwd <- rep(c(2, 1), each = 8)
    plotdata <- plotdata[, c("time", "NutrLum", "DLum", "RLum", "TransLum", "MdrLum", "MdtLum", "MrtLum", "MttLum",
                         "NutrWall", "DWall", "RWall", "TransWall", "MdrWall", "MdtWall", "MrtWall", "MttWall")]
  } else {
    mycol <- mycolother
    mylty <- myltyother
    mylwd <- rep(c(2, 1), each = 4)
    if(type == "Total") {
      maintitle <- "Pair model, totals"
      plotdata <- plotdata[, c("time", "NutrLum", "TotalDLum", "TotalRLum", "TotalTransLum",
                           "NutrWall", "TotalDWall", "TotalRWall", "TotalTransWall")]
    } else {
      maintitle <- "Bulk model"
    }
  }
  if(saveplot == TRUE) {
    filename <- paste0(DateTimeStamp, "output", gsub(" ", "", maintitle), ".png")
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      png(filename = filename)
    }
    matplot.deSolve(plotdata, main = maintitle,
                    sub = subtitle, ylim = myylim, log = if(yaxislog == 1) {"y"},
                    col = mycol, lty = mylty, lwd = mylwd,
                    legend = list(x = "bottomright"))
    grid()
    if(file.exists(filename) == FALSE) {
      dev.off()
    }
  } else {
    matplot.deSolve(plotdata, main = maintitle,
                    sub = subtitle, ylim = myylim, log = if(yaxislog == 1) {"y"},
                    col = mycol, lty = mylty, lwd = mylwd,
                    legend = list(x = "bottomright"))
    grid()
  }
}

# Parameterset 1: use wWall = wNutrWall to obtain equal recipient densities in
# the lumen and at the wall. Limit to wWall = wLum = wNutrLum to get equal
# nutrient concentrations in the lumen and at the wall. If those values are equal
# to the values used in the one-compartment model, then the densities and
# concentrations are also equal to those in the one-compartment model.
# This parameterset can be used to show influence of attachment and detachment
# rates that are different in the lumen compared to at the wall.
ScaleWallPerLum <- 1 # equal volumes in wall-compartment and lumen compartment
VLum <- 1
VWall <- 1
NILumSet <- 1.4
NIWallSet <- NILumSet
wLumSet <- -log(0.5)/24
wWallSet <- wLumSet
wNutrWallSet <- wLumSet
Ks <- 0.004
NutrConvSet <- 1.4e-7
bRSet <- 0.738
MigrLumWallSet <- 0.1
MigrWallLumSet <- MigrLumWallSet
DInitLumSet <- 1E3
DInitWallSet <- 0
cdSet <- 0.18
ctSet <- 0.09
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
kpWallSet <- c(0, 10^seq(from = -10, to = -8, by = 1))
knSet <- 10^seq(from = -1, to = 3, by = 0.1)
knWallSet <- 10^seq(from = -1, to = 3, length.out = 4)
gdSet <- 15
gtSet <- gdSet

# MyDataPar1 <- MyData

# Paramterset 1b, to get clearer plots
ScaleWallPerLum <- 1 # equal volumes in wall-compartment and lumen compartment
VLum <- 1
VWall <- 1
NILumSet <- 1.4
NIWallSet <- NILumSet
wLumSet <- -log(0.5)/24
wWallSet <- wLumSet
wNutrWallSet <- wLumSet
Ks <- 0.004
NutrConvSet <- 1.4e-7
bRSet <- 0.738
MigrLumWallSet <- 0.1
MigrWallLumSet <- MigrLumWallSet
DInitLumSet <- 1E3
DInitWallSet <- 0
cdSet <- 0.18
ctSet <- 0.09
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
kpWallSet <- kpSet
knSet <- 10^seq(from = -1, to = 3, length.out = 4)
knWallSet <- knSet
gdSet <- 15
gtSet <- gdSet

# Parameterset 2 to show effect of migration rates on biomass and stability of
# the plasmid-free equilibrium. Note that washout from the wall is excluded.
ScaleWallPerLum <- 1 # volume wall-compartiment / volume lumen compartiment
VLum <- 1
VWall <- 1
NILumSet <- 1.4
NIWallSet <- NILumSet
wLumSet <- -log(0.5)/24
wWallSet <- 0 # Cells do not washout from the wall
wNutrWallSet <- wLumSet
Ks <- 0.004
NutrConvSet <- 1.4e-7
bRSet <- 0.738
MigrLumWallSet <- c(0.025, 0.1, 0.4)
MigrWallLumSet <- MigrLumWallSet
DInitLumSet <- 1E3
DInitWallSet <- 0
cdSet <- 0.18
ctSet <- 0.09
kpSet <- 10^seq(from = -12, to = -8, by = 0.1)
kpWallSet <- 10^-12 
knSet <- 10^seq(from = -1, to = 3, by = 0.1)
knWallSet <- 10
gdSet <- 15
gtSet <- 15

# MyDataPar2 <- MyData

#### Main script ####
print(Sys.time())
CheckParms <- c(NILum = NILumSet, NIWall = NIWallSet, wLum = wLumSet,
                wWall = wWallSet, wNutrWallSet = wNutrWallSet,
                Ks = Ks, NutrConv = NutrConvSet, bR = bRSet,
                MigrLumWall = MigrLumWallSet, MigrWallLum = MigrWallLumSet,
                DInitLum = DInitLumSet, DInitWall = DInitWallSet, 
                cd = cdSet, ct = ctSet,
                kp = kpSet, kpWall = kpWallSet, kn = knSet, knWall = knWallSet,
                gd = gdSet, gt = gtSet)
warntext <- paste("Parameterset(s)",
                  paste(names(which(CheckParms <= 0)), collapse = ", "),
                  "contain(s) non-positive values.")
if(any(CheckParms <= 0)) warning(warntext)
if(any(c(cdSet, ctSet) >= 1)) warning("Costs should be larger than 0 and smaller than 1.")

TotalIterations <- length(NILumSet)*length(NIWallSet)*length(wLumSet)*
  length(wWallSet)*length(wNutrWallSet)*length(Ks)*length(NutrConvSet)*length(bRSet)*
  length(MigrLumWallSet)*length(MigrWallLumSet)*
  length(DInitLumSet)*length(DInitWallSet)*
  length(cdSet)*length(ctSet)*
  length(kpSet)*length(kpWallSet)*length(knSet)*length(knWallSet)*
  length(gdSet)*length(gtSet)
TotalIterations

## Get all parameter combinations to Determine plasmid-free equilibrium 
MyData <- expand_grid(NILum = NILumSet, NIWall = NIWallSet, wLum = wLumSet, wWall = wWallSet,
                      wNutrWall = wNutrWallSet, NutrConv = NutrConvSet, bR = bRSet,
                      MigrLumWall = MigrLumWallSet, MigrWallLum = MigrWallLumSet,
                      ScaleWallPerLum = ScaleWallPerLum)
dim(MyData)

# Multiply the equilibrium values of the one-compartment model with the
# appropriate volumes of the lumen or wall compartment to get milligram
# nutrients and number of recipient bacteria as guess for the plasmid-free
# equilibrium in the two-compartment model.
MyData <- mutate(MyData, 
                  NutrLumGuess = wLum*Ks*VLum / (bR - wLum),
                  RLumGuess = (NILum - NutrLumGuess)*VLum/NutrConv,
                  NutrWallGuess = wNutrWall*Ks*VWall / (bR - wWall),
                  RWallGuess = (NIWall - NutrWallGuess)*VWall/NutrConv
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
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")

## Approximate gdbulk and gtbulk in the lumen
VComp <- VLum
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
VComp <- VWall
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

# Approximate eigenvalues for pair-formation and bulk model
MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)
MyInfoEigVal <- t(apply(MyData, MARGIN = 1, FUN = CalcEigenvalues))
MyData <- cbind(MyData, MyInfoEigVal)

print("Eigenvalues estimated:")
print(Sys.time())

write.csv(MyData, file = paste0(DateTimeStamp, "outputnosimtwocomp.csv"),
          quote = FALSE, row.names = FALSE)

#### Create facet labels and labeller 'function' ####
labkpWall <- paste0("Attachment rate\nat the wall: ", signif(kpWallSet, 3))
names(labkpWall) <- kpWallSet
labknWall <- paste0("Detachment rate\nat the wall: ", signif(knWallSet, 3))
names(labknWall) <- knWallSet
labmigrlumwall <- paste0("Migration rate from\nlumen to wall: ", MigrLumWallSet)
names(labmigrlumwall) <- MigrLumWallSet
labmigrwalllum <- paste0("Migration rate from\nwall to lumen: ", MigrWallLumSet)
names(labmigrwalllum) <- MigrWallLumSet
if(length(knSet)==length(knWallSet)) {
  if(all(knSet == knWallSet)) {
    labkn <- paste0("Log10(detachment\nrates): ", signif(log10(knWallSet), 3))
    names(labkn) <- knSet
  }
}
if(exists("labkn")) {
  mylabeller <- labeller(kn = labkn, kpWall = labkpWall, knWall = labknWall,
                         MigrLumWall = labmigrlumwall,
                         MigrWallLum = labmigrwalllum, .default = label_both)  
} else {
  mylabeller <- labeller(kpWall = labkpWall, knWall = labknWall,
                         MigrLumWall = labmigrlumwall,
                         MigrWallLum = labmigrwalllum, .default = label_both)
}


#### Output parameterset 1 ####

# Show that biomass in lumen and in wall are equal
range(MyData$RLumInit/MyData$RWallInit)

# Show influence of kpWall and knWall on the stability of the plasmid-free
# equilibrium for parameterset 1 (Figure 5 in the article)
CreatePlot(filltitle = "Plasmid can invade",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Show influence of kpWall and knWall on the stability of the plasmid-free
# equilibrium in the bulk-model (not shown in article)
CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk-model)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Show if signs of bulk and pair-model are equal (Figure S4 in article)
CreatePlot(fillvar = "factor(SignDomEigVal == SignDomEigValBulk)",
           filltype = "manual",
           filltitle = "Dominant eigenvalues\nhave equal signs",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Recipient density in the lumen
CreatePlot(fillvar = "log10(RLumInit)", filltype = "continuous",
           filltitle = "Log10(Recipient density\nin the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Recipient density at the wall
CreatePlot(fillvar = "log10(RWallInit)", filltype = "continuous",
           filltitle = "Log10(Recipient density\nat the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Show that recipient density in lumen and in wall are equal
CreatePlot(fillvar = "RLumInit/RWallInit", filltype = "continuous",
           filltitle = "Recipient density\nin the lumen / at the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Bulk conjugation rates in the lumen
bulkrates <- c(MyData$gdbulkLum, MyData$gtbulkLum,
               MyData$gdbulkWall, MyData$gtbulkWall)
bulkrates <- bulkrates[which(bulkrates > 0)]
limitsbulkrates <- range(log10(bulkrates))

CreatePlot(fillvar = "log10(gdbulkLum)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate\nin the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(gtbulkLum)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant\nbulkrate in the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Bulk conjugation rates at the wall
CreatePlot(fillvar = "log10(gdbulkWall)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate\nat the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(gtbulkWall)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant\nbulkrate at the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Highest bulk-rate across lumen and wall
CreatePlot(fillvar = "log10(pmax(gdbulkLum, gdbulkWall))", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(max(Donor bulkrate\nin lumen and at wall))",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(pmax(gtbulkLum, gtbulkWall))", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(max(Transconjugant\nbulkrate in lumen and at the wall))",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Show percentage and counts of parameter combinations for which invasion is,
# or is not, possible
filteredDf <- NULL
for(i in kpWallSet) {
  for(j in knWallSet) {
    MyDataFiltered <- filter(MyData, kpWall == i, knWall == j)
    invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == 1))
    no_invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == -1))
    total_n <- invasion_n + no_invasion_n
    invasion_perc <- round(100*invasion_n/total_n, 0)
    no_invasion_perc <- round(100*no_invasion_n/total_n, 0)
    
    filteredDf <- rbind(filteredDf,
                        data.frame(kpWall = i, knWall = j,
                                   invasion_perc = invasion_perc,
                                   no_invasion_perc = no_invasion_perc,
                                   invasion_n = invasion_n,
                                   no_invasion_n = no_invasion_n, total_n = total_n))
  }
}
print(filteredDf)

# Alternative version (to be used with dataset 1b) to show percentage and counts
# of parameter combinations for which invasion is, or is not, possible
filteredDf <- NULL
for(i in knSet) {
  for(j in knWallSet) {
    MyDataFiltered <- filter(MyData, kn > i*0.999, kn < i*1.001,
                             knWall > j*0.999, knWall < j*1.001)
    invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == 1))
    no_invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == -1))
    total_n <- invasion_n + no_invasion_n
    invasion_perc <- round(100*invasion_n/total_n, 0)
    no_invasion_perc <- round(100*no_invasion_n/total_n, 0)
    
    filteredDf <- rbind(filteredDf,
                        data.frame(kn = i, knWall = j,
                                   invasion_perc = invasion_perc,
                                   no_invasion_perc = no_invasion_perc,
                                   invasion_n = invasion_n,
                                   no_invasion_n = no_invasion_n, total_n = total_n))
  }
}
print(filteredDf)


# To avoid another 4x4 plot, filter for facets where the detachment rate in the
# lumen equals the detachment rate at the wall.
FilteredData <- filter(MyData, kn == knWall)
ggplot(data = FilteredData,
       aes(x = log10(kp), y = log10(kpWall), fill = factor(SignDomEigVal))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  facet_wrap("kn", nrow = 2, labeller = mylabeller) +
  theme(legend.position = "bottom") +
  labs(x = "Log10(attachment rate in the lumen)",
       y = "Log10(attachment rate at the wall)", tag = NULL) +
  scale_fill_viridis_d("Plasmid can invade", labels = c("No", "Yes")) +
  geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1)
ggsave(paste0(DateTimeStamp, "InvasionTwoComp.png"))

ggplot(data = MyData,
       aes(x = log10(kp), y = log10(kpWall), fill = factor(SignDomEigVal))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  facet_grid(log10(kn) ~ log10(knWall), as.table = FALSE, labeller = mylabeller) +
  theme(legend.position = "bottom") +
  labs(x = "Log10(attachment rate in the lumen)",
       y = "Log10(attachment rate at the wall)", tag = NULL) +
  scale_fill_viridis_d("Plasmid can invade", labels = c("No", "Yes")) +
  geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1)
ggsave(paste0(DateTimeStamp, "InvasionTwoCompAlternative.png"))

ggplot(data = FilteredData,
       aes(x = log10(kp), y = log10(kpWall),
           fill = factor(SignDomEigVal==SignDomEigValBulk))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  facet_wrap("kn", nrow = 2, labeller = mylabeller) +
  theme(legend.position = "bottom") +
  labs(x = "Log10(attachment rate in the lumen)",
       y = "Log10(attachment rate at the wall)", tag = NULL) +
  scale_fill_manual("Dominant eigenvalues\nhave equal signs",
                    values = c("TRUE" = "darkgreen", "FALSE" = "red")) +
  geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1)
ggsave(paste0(DateTimeStamp, "InvasionTwoCompSameSign.png"))

ggplot(data = MyData,
       aes(x = log10(kp), y = log10(kpWall),
           fill = factor(SignDomEigVal==SignDomEigValBulk))) + 
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  facet_grid(log10(kn) ~ log10(knWall), as.table = FALSE, labeller = mylabeller) +
  theme(legend.position = "bottom") +
  labs(x = "Log10(attachment rate in the lumen)",
       y = "Log10(attachment rate at the wall)", tag = NULL) +
  scale_fill_manual("Dominant eigenvalues\nhave equal signs",
                    values = c("TRUE" = "darkgreen", "FALSE" = "red")) +
  geom_abline(intercept = 0, slope = 1, col = "white", size = 1.1)
ggsave(paste0(DateTimeStamp, "InvasionTwoCompAlternativeSameSign.png"))

# Separate plots for low attachment rates at the wall, to compare them to the
# output of the one-compartment model.
dim(MyData)
dim(MyData[which(MyData[, "kpWall"] == kpWallSet[1]), ])
DataWall1 <- MyData[which(MyData[, "kpWall"] == kpWallSet[1] & MyData[, "knWall"] == 10^-1), ]
DataWall2 <- MyData[which(MyData[, "kpWall"] == kpWallSet[1] & MyData[, "knWall"] == 10^1), ]
DataWall3 <- MyData[which(MyData[, "kpWall"] == kpWallSet[1] & MyData[, "knWall"] == 10^3), ]

CreatePlot(dataplot = DataWall1, filltitle = "Plasmid can invade",
           facetx = ".", facety = ".",
           filename = paste0(DateTimeStamp, "DataWall1.png"))
CreatePlot(dataplot = DataWall2, filltitle = "Plasmid can invade",
           facetx = ".", facety = ".",
           filename = paste0(DateTimeStamp, "DataWall2.png"))
CreatePlot(dataplot = DataWall3, filltitle = "Plasmid can invade",
           facetx = ".", facety = ".",
           filename = paste0(DateTimeStamp, "DataWall3.png"))

### Output to compare one-compartment model to two-compartment model

# Show that biomass in lumen and in wall are equal
range(MyData$RLumInit/MyData$RWallInit)

CreatePlot(filltitle = "Plasmid can invade",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk-model)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "factor(SignDomEigVal == SignDomEigValBulk)",
           filltype = "manual",
           filltitle = "Dominant eigenvalues\nhave equal signs",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

CreatePlot(fillvar = "log10(RLumInit)", filltype = "continuous",
           filltitle = "log10(Recipient density\nin the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(RWallInit)", filltype = "continuous",
           filltitle = "log10(Recipient density\nat the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "RLumInit/RWallInit", filltype = "continuous",
           filltitle = "Recipient density\nin the lumen / at the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

CreatePlot(fillvar = "log10(NutrLumInit)", filltype = "continuous",
           filltitle = "log10(Nutrient concentration\nin the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(NutrWallInit)", filltype = "continuous",
           filltitle = "log10(Nutrient concentration\nat the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "NutrLumInit/NutrWallInit", filltype = "continuous",
           filltitle = "Nutrient concentration\nin the lumen / at the wall",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Bulk conjugation rates in the lumen
bulkrates <- c(MyData$gdbulkLum, MyData$gtbulkLum,
               MyData$gdbulkWall, MyData$gtbulkWall)
bulkrates <- bulkrates[which(bulkrates > 0)]
limitsbulkrates <- range(log10(bulkrates))
CreatePlot(fillvar = "log10(gdbulkLum)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate\nin the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(gtbulkLum)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant\nbulkrate in the lumen)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Bulk conjugation rates at the wall
CreatePlot(fillvar = "log10(gdbulkWall)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate\nat the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(gtbulkWall)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant\nbulkrate at the wall)",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

# Highest bulk-rate across lumen and wall
CreatePlot(fillvar = "log10(pmax(gdbulkLum, gdbulkWall))", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(max(Donor bulkrate\nin lumen and at wall))",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)
CreatePlot(fillvar = "log10(pmax(gtbulkLum, gtbulkWall))", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(max(Transconjugant\nbulkrate in lumen and at the wall))",
           facetx = "kpWall", facety = "knWall", as.table = FALSE)

#### Output parameterset 2 ####

# Show the effect of migration rates on biomass at the wall (plot not shown).
CreatePlot(fillvar = "log10(RWallInit)", filltype = "continuous",
           filltitle = "Log10(Recipient\ndensity at the wall",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)

# Show the effect of migration rates on stability of the plasmid-free
# equilibrium (Figure 6).
CreatePlot(filltitle = "Plasmid can invade",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)

# Show the effect of migration rates on biomass in the lumen (plot not shown).
CreatePlot(fillvar = "log10(RLumInit)", filltype = "continuous",
           filltitle = "Log10(Recipient\ndensity in the lumen)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE,
           limits = range(log10(c(MyData$RLumInit, MyData$RWallInit))))

# Show the effect of migration rates on ratio of biomass at the wall to
# biomass in the lumen (plot not shown).
CreatePlot(fillvar = "log10(RWallInit/RLumInit)", filltype = "continuous",
           filltitle = "Log10(Recipient density at\nthe wall / in the lumen)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)

# Show the effect of migration rates on stability of the plasmid-free
# equilibrium for the bulk-model (plot not shown).
CreatePlot(fillvar = "factor(SignDomEigValBulk)",
           filltitle = "Plasmid can invade\n(bulk model)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)

# Are signs of the largest eigenvalues equal for bulk- and pair-formation model?
# (Figure S5 in article).
CreatePlot(fillvar = "factor(SignDomEigVal == SignDomEigValBulk)",
           filltype = "manual",
           filltitle = "Dominant eigenvalues\nhave equal signs",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)

# show effect of migration rates on the bulk-conjugation rates
limitsbulkrates <- range(log10(c(MyData$gdbulkLum, MyData$gtbulkLum,
                                 MyData$gdbulkWall, MyData$gtbulkWall)))
CreatePlot(fillvar = "log10(gdbulkLum)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor bulkrate\nin the lumen)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)
CreatePlot(fillvar = "log10(gtbulkLum)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant bulkrate\nin the lumen)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)
CreatePlot(fillvar = "log10(gdbulkWall)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Donor\nbulkrate at the wall)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)
CreatePlot(fillvar = "log10(gtbulkWall)", filltype = "continuous",
           limits = limitsbulkrates,
           filltitle = "Log10(Transconjugant\nbulkrate at the wall)",
           facetx = "MigrLumWall", facety = "MigrWallLum", as.table = FALSE)

MyDataTwoComp <- MyData
dim(MyDataTwoComp)
MyDataTwoComp1 <- filter(MyDataTwoComp, kpWall == kpWallSet[1], knWall == 10^-2)
MyDataTwoComp2 <- filter(MyDataTwoComp, kpWall == kpWallSet[1], knWall > 10^0.499 & knWall < 10^0.501)
MyDataTwoComp3 <- filter(MyDataTwoComp, kpWall == kpWallSet[1], knWall == 10^3)
dim(MyDataTwoComp1)

filteredDf <- NULL
for(i in MigrLumWallSet) {
  for(j in MigrWallLumSet) {
    MyDataFiltered <- filter(MyData, MigrLumWall == i, MigrWallLum == j)
    invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == 1))
    no_invasion_n <- length(which(MyDataFiltered[, "SignDomEigVal"] == -1))
    total_n <- invasion_n + no_invasion_n
    invasion_perc <- round(100*invasion_n/total_n, 0)
    no_invasion_perc <- round(100*no_invasion_n/total_n, 0)
    
    filteredDf <- rbind(filteredDf,
                        data.frame(MigrLumWall = i, MigrWallLum = j,
                                   invasion_perc = invasion_perc,
                                   no_invasion_perc = no_invasion_perc,
                                   invasion_n = invasion_n,
                                   no_invasion_n = no_invasion_n, total_n = total_n))
  }
}
print(filteredDf)


# See if bulk-rates differ for different detachment rates at the wall
CreatePlot(dataplot = MyDataTwoComp1, fillvar = "log10(gtbulkWall)", filltype = "continuous",
           limits = log10(range(MyDataTwoComp[, "gtbulkWall"])),
           filltitle = "Log10(Transconjugant\nbulkrate at the wall)",
           facetx = ".", facety = ".", save = FALSE)
CreatePlot(dataplot = MyDataTwoComp2, fillvar = "log10(gtbulkWall)", filltype = "continuous",
           limits = log10(range(MyDataTwoComp[, "gtbulkWall"])),
           filltitle = "Log10(Transconjugant\nbulkrate at the wall)",
           facetx = ".", facety = ".", save = FALSE)
CreatePlot(dataplot = MyDataTwoComp3, fillvar = "log10(gtbulkWall)", filltype = "continuous",
           limits = log10(range(MyDataTwoComp[, "gtbulkWall"])),
           filltitle = "Log10(Transconjugant\nbulkrate at the wall)",
           facetx = ".", facety = ".", save = FALSE)

CreatePlot(dataplot = filter(MyDataTwoComp, kpWall == 10^-12),
           filltitle = "Plasmid can invade",
           facetx = ".", facety = "knWall", save = FALSE)
CreatePlot(dataplot = MyDataTwoComp1,
           filltitle = "Plasmid can invade",
           facetx = ".", facety = ".", save = FALSE)
CreatePlot(dataplot = MyDataTwoComp2,
           filltitle = "Plasmid can invade",
           facetx = ".", facety = ".", save = FALSE)
CreatePlot(dataplot = MyDataTwoComp3,
           filltitle = "Plasmid can invade",
           facetx = ".", facety = ".", save = FALSE)
CreatePlot(dataplot = MyData,
           filltitle = "Plasmid can invade",
           facetx = ".", facety = ".", save = FALSE)

#### Final equilibria after invasion ####
# NOTE: running for kpLumSet = 10^seq(from = -12, to = -8, by = 0.1) and
# knLumSet = 10^seq(from = -2, to = 3, by = 0.1) leads to problems in the
# integration. Using by = 0.25 instead of by = 0.1 does work.

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
                              NutrLum = MyData[-IndexSimulation, "NutrLumInit"],
                              DLum = 0, RLum = MyData[-IndexSimulation, "RLumInit"],
                              TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0,
                              MttLum = 0,
                              NutrWall = MyData[-IndexSimulation, "NutrWallInit"],
                              DWall = 0, RWall = MyData[-IndexSimulation, "RWallInit"],
                              TransWall = 0, MdrWall = 0, MdtWall = 0, MrtWall = 0,
                              MttWall = 0)
  MyData <- rbind(cbind(MyData[IndexSimulation, ], OutputSimulationPairs),
                  cbind(MyData[-IndexSimulation, ], NoSimulationNeeded))
} else {
  MyData <- cbind(MyData, OutputSimulationPairs)
}
if(any(MyData$steady == 0)) warning("Steady-state has not always been reached")

print("Pair-formation model completed running:")
print(Sys.time())

write.csv(MyData, file = paste0(DateTimeStamp, "outputtwocompartmentpart1.csv"),
          quote = FALSE, row.names = FALSE)

IndexSimulationBulk <- which(MyData$SignDomEigValBulk != -1)
print(paste(length(IndexSimulationBulk), "simulations to run for the bulk model"))
InputSimulationBulk <- MyData[IndexSimulationBulk, ColumnsToSelect]
OutputSimulationBulk <- t(apply(X = InputSimulationBulk, MARGIN = 1, FUN = SimulationBulk))
colnames(OutputSimulationBulk) <- paste0(colnames(OutputSimulationBulk), "Bulk")

if(length(IndexSimulationBulk) < nrow(MyData)) {
  NoSimulationNeededBulk <- cbind(timeBulk = 0, steadyBulk = 1,
                                  NutrLumBulk = MyData[-IndexSimulation, "NutrLumInit"],
                                  DLumBulk = 0,
                                  RLumBulk = MyData[-IndexSimulation, "RLumInit"],
                                  TransLumBulk = 0,
                                  NutrWallBulk = MyData[-IndexSimulation, "NutrWallInit"],
                                  DWallBulk = 0,
                                  RWallBulk = MyData[-IndexSimulation, "RWallInit"],
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

##### Create plots over time #####
myltypairs <- c(lty = rep(c(3, 1, 2, 1, 1, 1, 1, 1), 2))
myltyother <- c(lty = rep(c(3, 1, 2, 1), 2))
mycolpairs <- rep(c("black", "purple", "green1", "red", "yellow", "hotpink", "blue", "cyan"), 2)
mycolother <- rep(c("black", "purple", "green1", "red"), 2)
myylim <- c(1E-6, 1E7)
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
verbose <- 0 # if verbose == 1, diagnositics on the simulations are printed
Mytmax <- c(500)
Mytstep <- c(0.1)
TheseRows <- c(1, nrow(MyData)) # Rows to use for simulations over time

ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))
Mydf <- MyData[TheseRows, ColumnsToSelect]
TotalIterations <- length(TheseRows)
print(TotalIterations)

# Times for which output of the simulation is wanted. Note that the used
# ode-solvers are variable-step methods, so the times in times are NOT the only
# times at which integration is performed. See help(diagnostics.deSolve()) and
# help(lsodar()) for details.
times <- seq(from = 0, to = Mytmax, by = Mytstep)

# To see the dynamics of the different populations
EqAfterInvasionPair <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime, type = "Pair"))
EqAfterInvasion <- cbind(Mydf, EqAfterInvasionPair)
for(i in TotalIterations) {dev.off()} # Not a clean way of implementing, but it works

# To compare total numbers of donors, recipients, and transconjugants in the
# output of the pair-formation model with the bulk-formation model
EqAfterInvasionTotal <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime, type = "Total"))
EqAfterInvasion <- cbind(Mydf, EqAfterInvasionTotal)
write.csv(EqAfterInvasion, file = paste0(DateTimeStamp, "outputtwocomprunovertime.csv"),
          quote = FALSE, row.names = FALSE)

#### The remainder of the script still has to be adjusted to the new model ####

########### Other plotting

limitsREq <- log10(c(min(c(MyData$RLumInit, MyData$RWallInit)), max(c(MyData$RLumInit, MyData$RWallInit))))
limitsNutrEq <- log10(c(min(c(MyData$NutrLumInit, MyData$NutrWallInit)), max(c(MyData$NutrLumInit, MyData$NutrWallInit))))

CreatePlot(fillvar = "log10(RLumInit)", limits = limitsREq, facetx = "NILum", facety = "NIWall")
CreatePlot2(fillvar = "log10(RLumInit)", filltype = "discrete", limits = c(-1, 1), facetx = "NILum", facety = "NIWall")

SetsWithMoreThanOneValue <- c("NILum", "NIWall", "wLum", "wWall", "NutrConv",
                              "bR", "MigrLumWall", "MigrWallLum")[c(length(NILumSet),
                             length(NIWallSet), length(wLumSet), length(wWallSet),
                             length(NutrConvSet), length(bRSet),
                             length(MigrLumWallSet), length(MigrWallLumSet)) > 1]

if(length(SetsWithMoreThanOneValue) < 5) {
  CreatePlot(fillvar = "RLumInit/RWallInit",
             xvar = SetsWithMoreThanOneValue[1],
             yvar = SetsWithMoreThanOneValue[2],
             facetx = ifelse(!is.na(SetsWithMoreThanOneValue[3]), SetsWithMoreThanOneValue[3], "."),
             facety = ifelse(!is.na(SetsWithMoreThanOneValue[4]), SetsWithMoreThanOneValue[4], "."))
}

if(length(SetsWithMoreThanOneValue) == 2) {
  CreatePlot(fillvar = "log10(RLumInit)", limits = limitsREq,
             xvar = SetsWithMoreThanOneValue[1],
             yvar = SetsWithMoreThanOneValue[2],
             facetx = SetsWithMoreThanOneValue[1], facety = SetsWithMoreThanOneValue[2])
}



if(length(SetsWithMoreThanOneValue) == 4) {
  CreatePlot(fillvar = "log10(RLumInit)", limits = limitsREq,
             xvar = "SetsWithMoreThanOneValue[1]",
             yvar = "SetsWithMoreThanOneValue[2]",
             facetx = "SetsWithMoreThanOneValue[3]", facety = "SetsWithMoreThanOneValue[4]")
}

CreatePlot2(fillvar = "log10(RLumInit)", limits = limitsREq, facetx = "SetsWithMoreThanOneValue", facety = "NIWall")

# CreatePlot(fillvar = "log10(RLumInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsREq)
# CreatePlot(fillvar = "log10(RWallInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsREq)
# CreatePlot(fillvar = "log10(NutrLumInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsNutrEq)
# CreatePlot(fillvar = "log10(NutrWallInit)", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall", limits = limitsNutrEq)
# 
# CreatePlot(fillvar = "RLumInit/RWallInit", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall")
# CreatePlot(fillvar = "NutrLumInit/NutrWallInit", xvar = "MigrLumWall", yvar = "MigrWallLum", facetx = "NILum", facety = "NIWall")



CreatePlot2(fillvar = "SignDomEigVal", gradient2 = 1, limits = c(-1, 1), facetx = "NILum", facety = "NIWall")

CreatePlot(fillvar = "SignDomEigVal", gradient2 = 1, limits = c(-1, 1), facetx = "MigrLumWall", facety = "MigrWallLum")
CreatePlot(fillvar = "SignDomEigVal", gradient2 = 1, limits = c(-1, 1), facetx = "NILum", facety = "NIWall")


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


CreatePlot(fillvar = "TotalRLum/TotalBioLum")
CreatePlot(fillvar = "TotalRWall/TotalBioWall")

CreatePlot(fillvar = "RLumBulk/TotalBioLumBulk")
CreatePlot(fillvar = "RWallBulk/TotalBioWallBulk")

CreatePlot(fillvar = "(RBulk/TotalBioBulk) / (TotalR/TotalBio)")
