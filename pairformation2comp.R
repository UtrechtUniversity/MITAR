#### Two-compartment pair-formation model of conjugation ####
# Copy-pasted the script of the one-compartment pairformation model.

# The upper part of the script should be able to run independently and applies to
# two-compartment pair-formation model of conjugation, as far as the change from
# the one-compartment to the two-compartment is implemented.

# The lower part of the script (after a line of #####) is the part of the
# script for the one-compartment model that is not yet changed to the
# two-compartment model.

### To do ###
# If I can prove that EqPlasmidfree1 is always non-positive I can remove the
# calculations and checks for it, and stick to using EqPlasmidfree2

#### Main script ####

#### Loading packages ####
library(deSolve) # Solve differential equations with results over time.
library(ggplot2) # For plotting data
library(RColorBrewer) # For better color schemes
library(rootSolve) # Integration, obtaining jacobian matrix and eigenvalues.
library(tidyr) # for 'expand.grid()' with dataframe as input

# Large set for testing
DInitSet <- c(1E3)
bRSet <- c(0.2, 0.8, 1.7)
NISet <- c(0.1, 1, 10)
NutrConv <- c(1e-8, 1e-6, 1e-4, 1e-2)
wSet <- c(0.01, 0.04, 0.10)
MigrLumWallSet <- c(0.01, 0.05, 0.1)
MigrWallLumSet <- c(0.01, 0.05, 0.1)
ScaleAreaPerVolSet <- c(0.5, 0.8, 1.5)
kpSet <- 10^seq(from = -10, to = -6, by = 0.5)
knSet <- 10^seq(from = -1, to = 3, by = 0.5)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(1, 15)
gtSet <- c(1, 15)

# Smaller set for testing
DInitSet <- c(1E3)
bRSet <- c(1.7)
NISet <- c(10)
NutrConv <- c(1e-6)
wSet <- c(0.04)
MigrLumWallSet <- c(0.05, 0.1)
MigrWallLumSet <- c(0.05, 0.1)
ScaleAreaPerVolSet <- c(0.8, 1.5)
kpSet <- 10^seq(from = -10, to = -6, by = 2)
knSet <- 10^seq(from = -1, to = 3, by = 2)
kpWallSet <- 10^seq(from = -10, to = -6, by = 2)
knWallSet <- 10^seq(from = -1, to = 3, by = 2)
cdSet <- c(0.01, 0.05)
ctSet <- c(0.01, 0.05)
gdSet <- c(15)
gtSet <- c(15)

CheckParms <- c(DInitSet, bRSet, NISet, NutrConv, MigrLumWallSet, MigrWallLumSet,
                ScaleAreaPerVolSet, wSet, kpSet, knSet, kpWallSet, knWallSet,
                cdSet, ctSet)
if(any(CheckParms <= 0)) warning("All parameters should have positive values.")
if(any(c(cdSet, ctSet) >= 1)) warning("Costs should be larger than 0 and smaller than 1.")

TotalIterations <- length(DInitSet)*length(bRSet)*length(NISet)*length(NutrConv)*
  length(wSet)*length(MigrLumWallSet)*length(MigrWallLumSet)*length(ScaleAreaPerVolSet)*
  length(kpSet)*length(knSet)*length(kpWallSet)*length(knWallSet)*length(cdSet)*
  length(ctSet)*length(gdSet)*length(gtSet)
TotalIterations

## Calculate plasmid-free equilibria for all parameter combinations
# Calculate the cell-containing plasmid-free equilibria (RLum*, RWall*, Nutr*),
# using the solution to 
# dRLum/dt = RLum*(Nutr*bR - w) - MigrLumWall*RLum + MigrWallLum*RWall*X == 0
# dRWall/dt = RWall*(Nutr*bR) + MigrLumWall*RLum/X - MigrWallLum*RWall == 0
# dNutr/dt = (NI - Nutr)*w - NutrConv*Nutr*bR*(RLum + X*RWall) == 0
# with RLum, RWall, Nutr > 0. The cell-free equilibrium is not considered but
# would be (RLum* = RWall* = 0, Nutr* = NI).
# I have used two functions for the solutions to make it easier to have Eq1 and
# Eq2 separate for further checks. Eq1 in calculations thus far always was
# non-positive (sometimes it was 0) and Eq2 usually (but not always) was 
# non-negative.
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

MyData <- expand_grid(bR = bRSet, NI = NISet, NutrConv = NutrConv, w = wSet,
                      MigrLumWall = MigrLumWallSet, MigrWallLum = MigrWallLumSet,
                      ScaleAreaPerVol = ScaleAreaPerVolSet)
dim(MyData)
Eqplasmidfree <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree1))
Eqplasmidfree2 <- t(apply(X = MyData, MARGIN = 1, FUN = CalcEqPlasmidfree2))

if(any(Eqplasmidfree[, "RLumEq"] <= 0 | Eqplasmidfree[, "RWallEq"] <= 0)) {
  IndexNegativeEq <- which(Eqplasmidfree[, "RLumEq"] <= 0 | Eqplasmidfree[, "RWallEq"] <= 0)
  ColnamesNegativeEq <- colnames(Eqplasmidfree)[unique(which(Eqplasmidfree <= 0, arr.ind = TRUE)[, 2])]
  warning("Equilibrium contains negative values in columns '", paste(ColnamesNegativeEq, collapse = "' and '"),
          "'.\nThe data will be included in the calculations anyway!
  This concerns the following rows of the dataframe: ", paste(IndexNegativeEq, collapse = ", "))
}

if(any(Eqplasmidfree2[, "RLumEq2"] > 0 & Eqplasmidfree2[, "RWallEq2"] > 0)) {
  IndexPositiveEq2 <- which(Eqplasmidfree2[, "RLumEq2"] > 0 & Eqplasmidfree2[, "RWallEq2"] > 0)
  ColnamesPositiveEq2 <- colnames(Eqplasmidfree2)[unique(which(Eqplasmidfree2 > 0, arr.ind = TRUE)[, 2])]
  warning("Equilibrium 2 contains positive values in columns '", paste(ColnamesPositiveEq2, collapse = "' and '"),
          "'.\nOnly equilibrium 1 will be used for calculations!
          This concerns the following rows of the dataframe: ", paste(ColnamesPositiveEq2, collapse = ", "))
}

head(Eqplasmidfree)
head(MyData)
MyData <- cbind(MyData, Eqplasmidfree)
head(MyData)

print("Plasmid-free equilibrium calculated:")
print(Sys.time())

## Add combinations with the parameters needed to approximate gdbulk and gtbulk
MyData <- expand_grid(MyData, kp = kpSet, kn = knSet, gd = gdSet, gt = gtSet,
                      DInit = DInitSet)
dim(MyData)

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
EstConjBulk <- function(MyData, RName) {
  with(as.list(MyData), {
    state <- c(D = MyData[["DInit"]], R = MyData[[RName]], Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0)
    parms <- MyData
    DataEstConjBulkDonor <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkDonor, parms = parms), 1)
    
    state <- c(R = MyData[[RName]], Trans = MyData[["DInit"]], Mrt = 0, Mtt = 0)
    DataEstConjBulkTrans <- tail(ode(t = timesEstConj, y = state,
                                     func = ModelEstConjBulkTrans, parms = parms), 1)
    
    DataEstConjBulk <- cbind(DataEstConjBulkDonor, DataEstConjBulkTrans)
    names(DataEstConjBulk) <- c(paste0("Donor", colnames(DataEstConjBulkDonor)),
                                paste0("Trans", colnames(DataEstConjBulkTrans)))
    return(DataEstConjBulk)
  })
}

timesEstConj <- seq(from = 0, to = 3, by = 0.1)
DataEstConjBulk <- t(apply(X = MyData, MARGIN = 1, FUN = EstConjBulk, RName = "RLumEq"))

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulk <- MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor)
gdbulk <- unname(gdbulk)
MyData <- cbind(MyData, gdbulkLum = gdbulk)

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] + DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] + DataEstConjBulk[, "TransMrt"]
gtbulk <- MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans)
gtbulk <- unname(gtbulk)
MyData <- cbind(MyData, gtbulk = gtbulk)

dim(MyData)
MyData <- expand_grid(MyData, kpWall = kpWallSet, knWall = knWallSet)
dim(MyData)

DataEstConjBulk <- t(apply(X = MyData, MARGIN = 1, FUN = EstConjBulk, RName = "RWallEq"))

TotalDEstConjBulkDonor <- DataEstConjBulk[, "DonorD"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMdt"]
TotalREstConjBulkDonor <- DataEstConjBulk[, "DonorR"] + DataEstConjBulk[, "DonorMdr"] + DataEstConjBulk[, "DonorMrt"]
gdbulk <- MyData[, "gd"] * DataEstConjBulk[, "DonorMdr"] / (TotalDEstConjBulkDonor * TotalREstConjBulkDonor)
gdbulk <- unname(gdbulk)
MyData <- cbind(MyData, gdbulkWall = gdbulk)

TotalTransEstConjBulkTrans <- DataEstConjBulk[, "TransTrans"] + DataEstConjBulk[, "TransMrt"] + 2*DataEstConjBulk[, "TransMtt"]
TotalREstConjBulkTrans <- DataEstConjBulk[, "TransR"] + DataEstConjBulk[, "TransMrt"]
gtbulk <- MyData[, "gt"] * DataEstConjBulk[, "TransMrt"] / (TotalTransEstConjBulkTrans * TotalREstConjBulkTrans)
gtbulk <- unname(gtbulk)
MyData <- cbind(MyData, gtbulkWall = gtbulk)

print("Bulk-conjugation rates estimated:")
print(Sys.time())

MyData <- expand_grid(MyData, cd = cdSet, ct = ctSet)


############################################################################################### 
##### Script for the one-compartment model, not yet adjusted to two-compartment situation #####
############################################################################################### 

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

# NOTE: hardcoded selection
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
                    name = "Plasmid-free equilibrium",
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

# If invasion is possible, run simulation to see how many bacteria of each
# population are present at equilibrium
IndexSimulation <- which(MyData$SignDomEigVal != -1)
print(paste(length(IndexSimulation), "simulations to run for the pair-formation model"))
ColumnsToSelect <- c(1:(which(names(MyData)=="Eigval1") - 1))
InputSimulationPairs <- MyData[IndexSimulation, ColumnsToSelect]
OutputSimulationPairs <- t(apply(X = InputSimulationPairs, MARGIN = 1,
                                 FUN = SimulationPairs))

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

IndexSimulationBulk <- which(MyData$SignDomEigValBulk != -1)
print(paste(length(IndexSimulation), "simulations to run for the bulk model"))
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

A <- MyData[which(MyData[, "TotalPlasmid"]/MyData[, "TotalBio"] > 1e-3 &
                    MyData[, "TotalPlasmid"]/MyData[, "TotalBio"] < 0.999), ]
dim(A)

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
BackupMyData <- MyData

###### To create plots over time

# Settings for simulations, plotting, and printing
mylty <- c(lty = c(3, 1, 2, 1, 1, 1, 1, 1))
# mycol <- c("black", "purple", "hotpink", "red", "yellow", "green1", "blue", "cyan")
mycol <- c("black", brewer.pal(7, "Set1"))
myylim <- c(1E-14, 1E7) # Defining the limits for the y-axis
yaxislog <- 1 # if yaxislog == 1, the y-axis is plotted on a logarithmic scale
plotoutput <- 1
extinctionthreshold <- 1E-10 # Population size is set to 0 if it is below the extinctionthreshold
verbose <- 0 # if verbose == 1, diagnositics on the simulations are printed and roots are indicated in the graphs
smallchange <- c(1E-5)
Mytmax <- c(1E4)
Mytstep <- c(10)

#### Create matrix to store data ####
TheseRows <- c(1, nrow(MyData))
Mydf <- MyData[TheseRows, ]
TotalIterations <- length(TheseRows)
print(TotalIterations)
CurrentIteration <- 0

# Times for which output of the simulation is wanted.
# Note that timesteps of 1 are used untill t = 100, and note furthermore that
# the used ode-solvers are variable-step methods, so the times in times
# are NOT the only times at which integration is performed. See
# help(diagnostics.deSolve()) and help(lsodar()) for details.
times <- c(0:100, seq(from = 100 + Mytstep, to = Mytmax, by = Mytstep))

# Define root-function and event-function
# If the root-argument is equal to 0, an event (as defined by the event-function)
# is triggered. See help(events) and help(lsodar) (both in the deSolve package)
# for background information and examples. 
# The first root becomes 0 if the sum of absolute rates of change is equal to the
# threshold smallchange, i.e., when equilibrium is nearly reached. This root is
# set as the terminal root, in order to terminate the integration if this occurs.
# The other roots are used to determine if any of the state variables are equal
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
rootfun <- function(t, state, parms) {
  c(sum(abs(unlist(ModelPairsNutr(t, state, parms)))) - smallchange, state - extinctionthreshold)
}

rootfunBulk <- function(t, stateBulk, parmsBulk) {
  c(sum(abs(unlist(ModelBulkNutr(t, stateBulk, parmsBulk)))) - smallchange, stateBulk - extinctionthreshold)
}

eventfun <- function(t, state, parms) {
  state[state < extinctionthreshold] <- 0
  return(state)
}

eventfunBulk <- function(t, stateBulk, parmsBulk) {
  stateBulk[stateBulk < extinctionthreshold] <- 0
  return(stateBulk)
}

RunOverTime <- function(data = Mydf, ...) {
  state <- c(Nutr = data[["NutrEq"]], D = data[["DInit"]], R = data[["REq"]],
             Trans = 0, Mdr = 0, Mdt = 0, Mrt = 0, Mtt = 0)
  parms <- data[ColumnsToSelect]
  out2 <- ode(t = times, y = state, func = ModelPairsNutr, parms = parms,
              rootfun = rootfun, events = list(func = eventfun, root = TRUE,
                                               terminalroot = 1))
  EqAfterInvasion <- tail(out2, 1)
  if(verbose == TRUE) {
    print(diagnostics(out2))
    print(attributes(out2))
  }
  PlotOverTime(data = out2, type = "Pair", saveplot = saveplots)
  
  stateBulk <- c(Nutr = data[["NutrEq"]], D = data[["DInit"]], R = data[["REq"]], Trans = 0)
  parmsBulk <- data[ColumnsToSelect]
  out2bulk <- ode(t = times, y = stateBulk, func = ModelBulkNutr, parms = parmsBulk,
                  rootfun = rootfunBulk, events = list(func = eventfunBulk, root = TRUE,
                                                       terminalroot = 1))
  EqAfterInvasionBulk <- tail(out2bulk, 1)
  if(verbose == TRUE) {
    print(diagnostics(out2bulk))
    print(attributes(out2bulk))
  }
  PlotOverTime(data = out2bulk, type = "Bulk", saveplot = saveplots)

  EqAfterInvasionTotal <- cbind(EqAfterInvasion, EqAfterInvasionBulk)
  names(EqAfterInvasionTotal) <- c(colnames(EqAfterInvasion),
                                   paste0(colnames(EqAfterInvasionBulk), "Bulk"))
  return(EqAfterInvasionTotal)
}

PlotOverTime <- function(data = out2, type = "Pair", saveplot = saveplots) {
  maintitle <- paste(type, "model")
  subtitle <- paste0("bR=", Mydf[["bR"]], " NI=", Mydf[["NI"]],
                     " log10kp=", log10(Mydf[["kp"]]), " log10kn=", log10(Mydf[["kn"]]),
                     " NutrConv=", Mydf[["NutrConv"]], " w=", Mydf[["w"]],
                     " cd=", Mydf[["cd"]], " ct=", Mydf[["ct"]])
  if(type == "Pair") {
    subtitle <- paste0(subtitle, " gd=", Mydf[["gd"]], " gt=", Mydf[["gt"]]) 
  } else {
    subtitle <- paste0(subtitle, " gdbulk=",  signif(Mydf[["gdbulk"]], digits = 4),
                       " gtbulk=", signif(Mydf[["gtbulk"]], digits = 4))
  }
  if(verbose == TRUE) {
    matplot.deSolve(data, main = paste(maintitle, "verbose"), sub = subtitle, ylim = myylim,
                    xlim = c(0, tail(attributes(data)$troot, 2)[1]),
                    log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                    legend = list(x = "topright"))
    grid()
    abline(h = extinctionthreshold)
    abline(v = attributes(data)$troot)
  }
  if(saveplot == TRUE) {
    filename <- paste0(DateTimeStamp, "output", gsub(" ", "", maintitle), ".png")
    if(file.exists(filename)) {
      warning("File already exists, not saved again!")
    } else {
      png(filename = filename)
    }
    matplot.deSolve(data, main = maintitle, sub = subtitle, ylim = myylim,
                    log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                    legend = list(x = "bottomright"))
    grid()
    if(saveplot == TRUE & file.exists(filename) == FALSE) {
      dev.off()
    }
  } else {
    matplot.deSolve(data, main = maintitle, sub = subtitle, ylim = myylim,
                    log = if(yaxislog == 1) {"y"}, col = mycol, lty = mylty, lwd = 2,
                    legend = list(x = "bottomright"))
    grid()
  }
}

EqAfterInvasionTotal <- t(apply(X = Mydf, MARGIN = 1, FUN = RunOverTime))


# png(filename = paste0(DateTimeStamp, "longrun", CurrentIteration + 1, ".png"))
# dev.off()

# MyData <- matrix(data = NA, nrow = TotalIterations, ncol = 61, byrow = TRUE) # To store output data


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
