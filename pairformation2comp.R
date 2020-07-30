#### Two-compartment pair-formation model of conjugation ####

#### To do ####
# Also see the 'To do' section in the pair-formation script for the one-compartment model

# Add CreatePlot() and SummaryPlot() function (see script for one-compartment
# model)

# Check if using MigrLumWall = MigrWallLum = 0 en DInitWall = 0 leads to same 
# results as the single-compartment model.

# If I can prove that EqPlasmidfree1 is always non-positive I can remove the
# calculations and checks for it, and stick to using EqPlasmidfree2

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

#### Plotting options ####
saveplots <- 0
tmaxsteady <- 1e8
timesEstConj <- seq(from = 0, to = 3, by = 0.1)
MyColorBrew <- rev(brewer.pal(11, "Spectral")) # examples: display.brewer.all()
# mycol <- c("black", brewer.pal(7, "Set1"))

#### Functions ####

# Calculate the cell-containing plasmid-free equilibria (RLum*, RWall*, Nutr*),
# using the solution to 
# dRLum/dt = RLum*(Nutr*bR - w) - MigrLumWall*RLum + MigrWallLum*RWall*X == 0
# dRWall/dt = RWall*(Nutr*bR) + MigrLumWall*RLum/X - MigrWallLum*RWall == 0
# dNutr/dt = (NI - Nutr)*w - NutrConv*Nutr*bR*(RLum + X*RWall) == 0
# with RLum, RWall, Nutr > 0. The cell-free equilibrium is not considered but
# would be (RLum* = RWall* = 0, Nutr* = NI).
# I have used two functions for the solutions to make it easier to have Eq1 and
# Eq2 separate for further checks. Eq1 in calculations thus far usually (but not
# always) was non-negative and Eq2 always was non-positive (sometimes it was 0).
CalcEqPlasmidfree1 <- function(MyData) {
  with(as.list(MyData), {
    NutrEq <- (MigrLumWall + MigrWallLum + w -
                 (MigrLumWall^2 + 2*MigrLumWall*MigrWallLum + 2*MigrLumWall*w +
                    MigrWallLum^2 - 2*MigrWallLum*w + w^2)^(1/2))/(2*bR)
    RLumEq <- NI/NutrConv - (MigrLumWall + MigrWallLum + w -
                               (MigrLumWall^2 + 2*MigrLumWall*MigrWallLum +
                                  2*MigrLumWall*w + MigrWallLum^2 -
                                  2*MigrWallLum*w + w^2)^(1/2))/(2*NutrConv*bR)
    RWallEq <- (MigrLumWall*NI*bR -
                  MigrWallLum*w + NI*bR*w)/(MigrWallLum*NutrConv*ScaleAreaPerVol*bR) +
      ((MigrWallLum - NI*bR)*(MigrLumWall + MigrWallLum + w -
                                (MigrLumWall^2 + 2*MigrLumWall*MigrWallLum + 2*MigrLumWall*w +
                                   MigrWallLum^2 - 2*MigrWallLum*w + w^2)^(1/2)))/(2*MigrWallLum*NutrConv*ScaleAreaPerVol*bR)
    Eq <- c(NutrEq = NutrEq, RLumEq = RLumEq, RWallEq = RWallEq)
    return(Eq)
  })
}

CalcEqPlasmidfree2 <- function(MyData) {
  with(as.list(MyData), {
    NutrEq <- (MigrLumWall + MigrWallLum + w +
                 (MigrLumWall^2 + 2*MigrLumWall*MigrWallLum + 2*MigrLumWall*w +
                    MigrWallLum^2 - 2*MigrWallLum*w + w^2)^(1/2))/(2*bR)
    RLumEq <- NI/NutrConv - (MigrLumWall + MigrWallLum + w +
                               (MigrLumWall^2 + 2*MigrLumWall*MigrWallLum +
                                  2*MigrLumWall*w + MigrWallLum^2 -
                                  2*MigrWallLum*w + w^2)^(1/2))/(2*NutrConv*bR)
    RWallEq <- (MigrLumWall*NI*bR -
                  MigrWallLum*w + NI*bR*w)/(MigrWallLum*NutrConv*ScaleAreaPerVol*bR) +
      ((MigrWallLum - NI*bR)*(MigrLumWall + MigrWallLum + w +
                                (MigrLumWall^2 + 2*MigrLumWall*MigrWallLum + 2*MigrLumWall*w +
                                   MigrWallLum^2 - 2*MigrWallLum*w + w^2)^(1/2)))/(2*MigrWallLum*NutrConv*ScaleAreaPerVol*bR)
    Eq2 <- c(NutrEq2 = NutrEq, RLumEq2 = RLumEq, RWallEq2 = RWallEq)
    return(Eq2)
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
      warning("DInitLum == 0, using DInitWall for approximation of bulkrates in the lumen instead!")
      state <- c(D = MyData[["DInitWall"]], R = MyData[["RLumEq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    } else {
      state <- c(D = MyData[["DInitLum"]], R = MyData[["RLumEq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    }
    
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    if(MyData[["DInitLum"]] == 0) {
      state <- c(R = MyData[["RLumEq"]], Trans = MyData[["DInitWall"]], Mrt = 0, Mtt = 0)
    } else {
      state <- c(R = MyData[["RLumEq"]], Trans = MyData[["DInitLum"]], Mrt = 0, Mtt = 0)
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
      warning("DInitWall == 0, using DInitLum for approximation of bulkrates at the wall instead!")
      state <- c(D = MyData[["DInitLum"]], R = MyData[["RWallEq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    } else {
      state <- c(D = MyData[["DInitWall"]], R = MyData[["RWallEq"]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    }
    
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    if(MyData[["DInitWall"]] == 0) {
      state <- c(R = MyData[["RWallEq"]], Trans = MyData[["DInitLum"]], Mrt = 0, Mtt = 0)
    } else {
      state <- c(R = MyData[["RWallEq"]], Trans = MyData[["DInitWall"]], Mrt = 0, Mtt = 0)
    }
    
    DataEstConjBulkTrans <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)
    
    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
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
# Imran (2005).
ModelPairsNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w -
      NutrConv*Nutr*(
        (1 - cd)*bR*(DLum + MdrLum + MdtLum + ScaleAreaPerVol*(DWall + MdrWall + MdtWall)) +
          bR*(RLum + MdrLum + MrtLum + ScaleAreaPerVol*(RWall + MdrWall + MrtWall)) +
          (1 - ct)*bR*(TransLum + MdtLum + MrtLum + 2*MttLum +
                         ScaleAreaPerVol*(TransWall + MdtWall + MrtWall + 2*MttWall))
      )
    
    dDLum <- (1 - cd)*bR*Nutr*(DLum + MdrLum + MdtLum) - kp*DLum*RLum +
      kn*(MdrLum + MdtLum) - (w + MigrLumWall)*DLum + MigrWallLum*DWall*ScaleAreaPerVol
    dRLum <- bR*Nutr*(RLum + MdrLum + MrtLum) - kp*RLum*(DLum + TransLum) +
      kn*(MdrLum + MrtLum) - (w + MigrLumWall)*RLum + MigrWallLum*RWall*ScaleAreaPerVol
    dTransLum <- (1 - ct)*bR*Nutr*(TransLum + MdtLum + MrtLum + 2*MttLum) -
      kp*RLum*TransLum + kn*(MdtLum + MrtLum + 2*MttLum) -
      (w + MigrLumWall)*TransLum + MigrWallLum*TransWall*ScaleAreaPerVol
    dMdrLum <- kp*DLum*RLum - kn*MdrLum - gd*MdrLum - (w + MigrLumWall)*MdrLum +
      MigrWallLum*MdrWall*ScaleAreaPerVol
    dMdtLum <- gd*MdrLum - kn*MdtLum - (w + MigrLumWall)*MdtLum +
      MigrWallLum*MdtWall*ScaleAreaPerVol
    dMrtLum <- kp*RLum*TransLum - kn*MrtLum - gt*MrtLum -
      (w + MigrLumWall)*MrtLum + MigrWallLum*MrtWall*ScaleAreaPerVol
    dMttLum <- gt*MrtLum - kn*MttLum -
      (w + MigrLumWall)*MttLum + MigrWallLum*MttWall*ScaleAreaPerVol
    
    dDWall <- (1 - cd)*bR*Nutr*(DWall + MdrWall + MdtWall) - kpWall*DWall*RWall +
      knWall*(MdrWall + MdtWall) + MigrLumWall*DLum/ScaleAreaPerVol - MigrWallLum*DWall
    dRWall <- bR*Nutr*(RWall + MdrWall + MrtWall) - kpWall*RWall*(DWall + TransWall) +
      knWall*(MdrWall + MrtWall) + MigrLumWall*RLum/ScaleAreaPerVol - MigrWallLum*RWall
    dTransWall <- (1 - ct)*bR*Nutr*(TransWall + MdtWall + MrtWall + 2*MttWall) -
      kpWall*RWall*TransWall + knWall*(MdtWall + MrtWall + 2*MttWall) +
      MigrLumWall*TransLum/ScaleAreaPerVol - MigrWallLum*TransWall
    dMdrWall <- kpWall*DWall*RWall - knWall*MdrWall - gd*MdrWall +
      MigrLumWall*MdrLum/ScaleAreaPerVol - MigrWallLum*MdrWall
    dMdtWall <- gd*MdrWall - knWall*MdtWall + MigrLumWall*MdtLum/ScaleAreaPerVol - 
      MigrWallLum*MdtWall
    dMrtWall <- kpWall*RWall*TransWall - knWall*MrtWall - gt*MrtWall + 
      MigrLumWall*MrtLum/ScaleAreaPerVol - MigrWallLum*MrtWall
    dMttWall <- gt*MrtWall - knWall*MttWall + 
      MigrLumWall*MttLum/ScaleAreaPerVol - MigrWallLum*MttWall
    
    return(list(c(dNutr, dDLum, dRLum, dTransLum, dMdrLum, dMdtLum, dMrtLum, dMttLum,
                  dDWall, dRWall, dTransWall, dMdrWall, dMdtWall, dMrtWall, dMttWall)))
  })
}

# Bulk-conjugation model
ModelBulkNutr <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dNutr <- (NI - Nutr)*w -
      NutrConv*Nutr*(
        (1 - cd)*bR*(DLum + ScaleAreaPerVol*DWall) +
          bR*(RLum + ScaleAreaPerVol*RWall) +
          (1 - ct)*bR*(TransLum + ScaleAreaPerVol*TransWall) 
      )
    
    dDLum <- (1 - cd)*bR*Nutr*DLum - (w + MigrLumWall)*DLum +
      MigrWallLum*DWall*ScaleAreaPerVol
    dRLum <- bR*Nutr*RLum - gdbulkLum*DLum*RLum - gtbulkLum*TransLum*RLum - 
      (w + MigrLumWall)*RLum + MigrWallLum*RWall*ScaleAreaPerVol
    dTransLum <- (1 - ct)*bR*Nutr*TransLum + gdbulkLum*DLum*RLum + 
      gtbulkLum*TransLum*RLum - (w + MigrLumWall)*TransLum + MigrWallLum*TransWall*ScaleAreaPerVol
    
    dDWall <- (1 - cd)*bR*Nutr*DWall - MigrWallLum*DWall + MigrLumWall*DLum/ScaleAreaPerVol
    dRWall <- bR*Nutr*RWall - gdbulkWall*DWall*RWall - gtbulkWall*TransWall*RWall -
      MigrWallLum*RWall + MigrLumWall*RLum/ScaleAreaPerVol
    dTransWall <- (1 - ct)*bR*Nutr*TransWall + gdbulkWall*DWall*RWall +
      gtbulkWall*TransWall*RWall - MigrWallLum*TransWall + MigrLumWall*TransLum/ScaleAreaPerVol
    
    return(list(c(dNutr, dDLum, dRLum, dTransLum, dDWall, dRWall, dTransWall)))
  })
}

# Numerically estimate the Jacobian matrix of the plasmid-free equilibrium of
# the models, then calculate (or approximate?) the eigenvalues of this matrix.
# The maximum real part of the eigenvalues is used to determine stability.
CalcEigenvalues <- function(MyData) {
  parms <- MyData
  EqFull <- c(Nutr = MyData[["NutrEq"]], DLum = 0, RLum = MyData[["RLumEq"]],
              TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0, MttLum = 0,
              DWall = 0, RWall = MyData[["RWallEq"]],
              TransWall = 0, MdrWall = 0, MdtWall = 0, MrtWall = 0, MttWall = 0)
  EigVal <- eigen(x = jacobian.full(y = EqFull, func = ModelPairsNutr, parms = parms),
                  symmetric = FALSE, only.values = TRUE)$values
  ComplexEigVal <- is.complex(EigVal) 
  EigVal <- Re(EigVal)
  names(EigVal) <- paste0("Eigval", 1:length(EigVal))
  DomEigVal <- max(EigVal)
  SignDomEigVal <- sign(DomEigVal)
  SignEigValEqual <- identical(rep(SignDomEigVal, length(EigVal)), sign(Re(EigVal)))
  
  EqFullBulk <- c(Nutr = MyData[["NutrEq"]], DLum = 0, RLum = MyData[["RLumEq"]], TransLum = 0, 
                  DWall = 0, RWall = MyData[["RWallEq"]], TransWall = 0)
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
# as state. Note that stol is based on the average of absolute rates of change,
# not the sum.
SimulationPairs <- function(InputSimulationPairs) {
  parms <- InputSimulationPairs
  state <- c(Nutr = parms[["NutrEq"]], DLum = parms[["DInitLum"]],
             RLum = parms[["RLumEq"]], TransLum = 0, MdrLum = 0, MdtLum = 0,
             MrtLum = 0, MttLum = 0, DWall = parms[["DInitWall"]],
             RWall = parms[["RWallEq"]], TransWall = 0, MdrWall = 0, MdtWall = 0,
             MrtWall = 0, MttWall = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelPairsNutr,
                   parms = parms, stol = 1.25e-6)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

SimulationBulk <- function(InputSimulationBulk) {
  parms <- InputSimulationBulk
  state <- c(Nutr = parms[["NutrEq"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumEq"]], TransLum = 0, 
             DWall = parms[["DInitWall"]], RWall = parms[["RWallEq"]], TransWall = 0)
  out <- runsteady(y = state, time = c(0, tmaxsteady), func = ModelBulkNutr, parms = parms, stol = 2.5e-6)
  EqAfterInvDonor <- c(time = attr(out, "time"), steady = attr(out, "steady"), out$y)
  return(EqAfterInvDonor)
}

RunOverTime <- function(parms = Mydf, verbose = FALSE, type = "Pair", ...) {
  state <- c(Nutr = parms[["NutrEq"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumEq"]],
             TransLum = 0, MdrLum = 0, MdtLum = 0, MrtLum = 0, MttLum = 0,
             DWall = parms[["DInitWall"]], RWall = parms[["RWallEq"]],
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
  stateBulk <- c(Nutr = parms[["NutrEq"]], DLum = parms[["DInitLum"]], RLum = parms[["RLumEq"]], TransLum = 0,
                 DWall = parms[["DInitWall"]], RWall = parms[["RWallEq"]], TransWall = 0)
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

#### Parameter values ####
# MigrLumWallSet: rate of migration from the lumen to the wall
# MigrWallLumSet: rate of migration from the wall to the lumen
# See script of the one-compartment model for explanation of the other parameter values

# Large set for testing
DInitLumSet <- c(1E3)
DInitWallSet <- c(1E3)
bRSet <- c(0.2, 0.8, 1.7)
NISet <- c(0.1, 1, 10)
NutrConvSet <- c(1e-8, 1e-6, 1e-4, 1e-2)
wSet <- c(0.01, 0.04, 0.10)
MigrLumWallSet <- c(0.01, 0.05, 0.1)
MigrWallLumSet <- c(0.01, 0.05, 0.1)
ScaleAreaPerVolSet <- c(0.5, 0.8, 1.5)
kpSet <- 10^seq(from = -10, to = -6, by = 0.5)
knSet <- 10^seq(from = -1, to = 3, by = 0.5)
kpWallSet <- 10^seq(from = -10, to = -6, by = 2)
knWallSet <- 10^seq(from = -1, to = 3, by = 2)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(1, 15)
gtSet <- c(1, 15)

# Smaller set for testing
DInitLumSet <- c(1E3)
DInitWallSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
NutrConvSet <- c(1e-6)
wSet <- c(0.04)
MigrLumWallSet <- c(0.05)
MigrWallLumSet <- c(0.1)
ScaleAreaPerVolSet <- c(0.8)
kpSet <- 10^seq(from = -10, to = -6, by = 2)
knSet <- 10^seq(from = -1, to = 1, by = 2)
# kpSet <- 10^-8
# knSet <- 10
# kpWallSet <- 10^seq(from = -10, to = -6, by = 2)
# knWallSet <- 10^seq(from = -1, to = 3, by = 2)
kpWallSet <- 10^-8
knWallSet <- 1
cdSet <- c(0.05)
ctSet <- c(0.05)
gdSet <- c(15)
gtSet <- c(15)


# Testing influence of different REq on approximations of bulk conjugationrates
# Without migration from lumen to wall or vice versa
DInitLumSet <- c(1E3)
DInitWallSet <- c(0)
bRSet <- c(1.7)
NISet <- c(1, 10, 100)
NutrConvSet <- c(1e-6)
wSet <- c(0.04)
MigrLumWallSet <- c(0.05)
MigrWallLumSet <- c(0.1)
ScaleAreaPerVolSet <- c(1)
kpSet <- 10^-8
knSet <- 1
kpWallSet <- 10^-8
knWallSet <- 1
cdSet <- c(0.05)
ctSet <- c(0.05)
gdSet <- c(15)
gtSet <- c(15)

## Testset for code testing
DInitLumSet <- c(1E3)
DInitWallSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(1, 10, 100)
NutrConvSet <- c(1e-6)
wSet <- c(0.04)
MigrLumWallSet <- c(0.05)
MigrWallLumSet <- c(0.1)
ScaleAreaPerVolSet <- c(1)
kpSet <- 10^-8
knSet <- 1
kpWallSet <- c(10^-8, 10^-6)
knWallSet <- c(1, 10)
cdSet <- c(0.05)
ctSet <- c(0.05)
gdSet <- c(15)
gtSet <- c(15)


### Why do DLum and DWall behave differently from each other in the first ~200h hours
### in the pair-formation model (the drop in densities is not seen in bulk model)
DInitLumSet <- c(1E3)
DInitWallSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
NutrConvSet <- c(1e-6)
wSet <- c(0.04)
MigrLumWallSet <- c(0.05)
MigrWallLumSet <- c(0.05)
ScaleAreaPerVolSet <- c(1.2)
kpSet <- 10^-10
knSet <- c(0.1, 10)
kpWallSet <- c(10^-8, 10^-6)
knWallSet <- 0.1
cdSet <- c(0.05)
ctSet <- c(0.05)
gdSet <- c(15)
gtSet <- c(15)

#### Main script ####
CheckParms <- c(DInitLumSet = DInitLumSet, DInitWallSet = DInitWallSet, bRSet = bRSet,
                NISet = NISet, NutrConvSet = NutrConvSet, MigrLumWallSet = MigrLumWallSet,
                MigrWallLumSet = MigrWallLumSet, ScaleAreaPerVolSet = ScaleAreaPerVolSet,
                wSet = wSet, kpSet = kpSet, knSet = knSet, kpWallSet = kpWallSet,
                knWallSet = knWallSet, cdSet = cdSet, ctSet = ctSet)
if(any(CheckParms <= 0)) warning(cat("Parameter(s)", names(which(CheckParms <= 0)), "contain(s) non-positive values."))
if(any(c(cdSet, ctSet) >= 1)) warning("Costs should be larger than 0 and smaller than 1.")

TotalIterations <- length(DInitLumSet)*length(DInitWallSet)*length(bRSet)*length(NISet)*length(NutrConvSet)*
  length(wSet)*length(MigrLumWallSet)*length(MigrWallLumSet)*length(ScaleAreaPerVolSet)*
  length(kpSet)*length(knSet)*length(kpWallSet)*length(knWallSet)*length(cdSet)*
  length(ctSet)*length(gdSet)*length(gtSet)
TotalIterations

## Calculate plasmid-free equilibrium for all parameter combinations
MyData <- expand_grid(bR = bRSet, NI = NISet, NutrConv = NutrConvSet, w = wSet,
                      MigrLumWall = MigrLumWallSet, MigrWallLum = MigrWallLumSet,
                      ScaleAreaPerVol = ScaleAreaPerVolSet)
dim(MyData)
Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree1))
Eqplasmidfree2 <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree2))

if(any(Eqplasmidfree <= 0)) {
  RowsNegativeEq <- unique(which(Eqplasmidfree <= 0, arr.ind = TRUE)[, 1])
  ColnamesNegativeEq <- colnames(Eqplasmidfree)[unique(which(Eqplasmidfree <= 0, arr.ind = TRUE)[, 2])]
  warning("Equilibrium contains non-positive values in columns '", paste(ColnamesNegativeEq, collapse = "' and '"),
          "'.\nThe data will be included in the calculations anyway!
  This concerns the following rows of the dataframe: ", paste(RowsNegativeEq, collapse = ", "))
}

if(any(Eqplasmidfree2[, "RLumEq2"] > 0 & Eqplasmidfree2[, "RWallEq2"] > 0)) {
  RowsPositiveEq2 <- which(Eqplasmidfree2[, "RLumEq2"] > 0 & Eqplasmidfree2[, "RWallEq2"] > 0)
  ColnamesPositiveEq2 <- colnames(Eqplasmidfree2)[unique(which(Eqplasmidfree2 > 0, arr.ind = TRUE)[, 2])]
  warning("Equilibrium 2 contains positive values in columns '", paste(ColnamesPositiveEq2, collapse = "' and '"),
          "'.\nOnly equilibrium 1 will be used for calculations!
          This concerns the following rows of the dataframe: ", paste(RowsPositiveEq2, collapse = ", "))
}

MyData <- cbind(MyData, Eqplasmidfree)
print("Plasmid-free equilibrium calculated:")
print(Sys.time())

## Approximate gdbulk and gtbulk in the lumen
MyData <- expand_grid(MyData, gd = gdSet, gt = gtSet, DInitLum = DInitLumSet, DInitWall = DInitWallSet, kp = kpSet, kn = knSet)
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

# Replace columns kn and kn with the values in the columns kpWall and knWall in new dataframe
# to be used to estimate bulk conjugation rates at the wall
MyDataWall <- cbind(MyData[, 1:(which(names(MyData)=="kp") - 1)],
                    kp = unname(MyData[, "kpWall"]), kn = unname(MyData[, "knWall"]))
DataEstConjBulk <- t(apply(X = MyDataWall, MARGIN = 1, FUN = EstConjBulkWall))

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulkWall <- unname(MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor))

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] + DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] + DataEstConjBulk[, "TransMrt"]
gtbulkWall <- unname(MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans))

MyData <- cbind(MyData, gdbulkWall = gdbulkWall, gtbulkWall = gtbulkWall)

print("Bulk-conjugation rates estimated:")
print(Sys.time())

# calculate (or approximate?) eigenvalues
MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)

MyInfoEigVal <- t(apply(MyData, MARGIN = 1, FUN = CalcEigenvalues))
MyData <- cbind(MyData, MyInfoEigVal)

print("Eigenvalues estimated:")
print(Sys.time())

DateTimeStamp <- format(Sys.time(), format = "%Y_%B_%d_%H_%M_%S")
write.csv(MyData, file = paste0(DateTimeStamp, "outputnosimulationtwocompartment.csv"),
          quote = FALSE, row.names = FALSE)

# If invasion is possible, run simulation to see how many bacteria of each
# population are present at equilibrium

IndexSimulation <- which(MyData$SignDomEigVal != -1)
print(paste(length(IndexSimulation), "simulations to run for the pair-formation model"))
ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))
InputSimulationPairs <- MyData[IndexSimulation, ColumnsToSelect]
OutputSimulationPairs <- t(apply(X = InputSimulationPairs, MARGIN = 1,
                                 FUN = SimulationPairs))

if(length(IndexSimulation) < nrow(MyData)) {
  NoSimulationNeeded <- cbind(time = 0, steady = 1, Nutr = MyData[-IndexSimulation, "NutrEq"],
                              DLum = 0, RLum = MyData[-IndexSimulation, "RLumEq"],
                              TransLum = 0,  MdrLum = 0, MdtLum = 0, MrtLum = 0,
                              MttLum = 0, DWall = 0,
                              RWall = MyData[-IndexSimulation, "RWallEq"],
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
print(paste(length(IndexSimulation), "simulations to run for the bulk model"))
InputSimulationBulk <- MyData[IndexSimulationBulk, ColumnsToSelect]
OutputSimulationBulk <- t(apply(X = InputSimulationBulk, MARGIN = 1, FUN = SimulationBulk))
colnames(OutputSimulationBulk) <- paste0(colnames(OutputSimulationBulk), "Bulk")

if(length(IndexSimulationBulk) < nrow(MyData)) {
  NoSimulationNeededBulk <- cbind(timeBulk = 0, steadyBulk = 1,
                                  NutrBulk = MyData[-IndexSimulation, "NutrEq"],
                                  DLumBulk = 0,
                                  RLumBulk = MyData[-IndexSimulation, "RLumEq"],
                                  TransLumBulk = 0, DWallBulk = 0,
                                  RWallBulk = MyData[-IndexSimulation, "RWallEq"],
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
