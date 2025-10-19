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


#### References ####
# Roberts MG, Heesterbeek JAP. 2021. Infection dynamics in ecosystems: on the
# interaction between red and grey squirrels, pox virus, pine martens and trees.
# Journal of the Royal Society, Interface 18(183):20210551.

# Wickham H, Grolemund G. 2017. R for data science: import, tidy, transform,
# visualize, and model data. Online version: https://r4ds.had.co.nz/index.html


#### Loading required packages ####
library(dplyr)     # across(), full_join(), group_by(), near(), summarise()
library(ggplot2)   # to display data and results
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()
# On the pipe operator (%>%), see ?'%>%' and Ch. 18 'Pipes' in Wickham 2017.


#### Load required functions ####
source("./ms_funcs.R")


#### Settings and defining parameter space ####
# Note: to simulate that each species belongs to a different class, an
# additional taxmattype has to be added. Then only the interspecies conjugation
# rates should be reduced, as the intraspecies conjugation rates of all species
# are still those given in conjrateset. This unchanged intraspecies conjugation
# rate makes it different from reducing conjrateset 1000-fold and using
# 'taxmatsame' as taxmattype.

## Basis parameter set to create bifurcation-like plots showing the border of
# epidemiological stability in the conjugation rate/cost space
saveplots <- TRUE
niterintmat <- 1
smallstate <- NA
finalsmallstate <- NA
smallchange <- NA
totalabun <- 1e11
nspeciesset <- c(2, 16)
maxnspecies <- max(nspeciesset)
abunmodelset <- c("dompreempt")
costtype <- "absolute"
costmark <- NULL # Plot dotted vertical lines at indicated values if not NULL
conjrate_base <- 1e-12
seqconjrate <- 10^seq(from = -13.0, to = -11.5, by = 0.005)
# If taxmattype is "SameSpecies", the conjugation rate is the same for all
# populations. If taxmattype is "OtherClass", the interspecies conjugation rate
# to and from the initially plasmid-bearing population is reduced a 1000-fold.
taxmattypeset <- c("SameSpecies", "OtherClass")
# Some plasmid-free bacteria are added to simulate perturbation by plasmid-free
# bacteria, and some plasmid-free bacteria are replaced with plasmid-bearing
# bacteria to simulate perturbation by plasmid-bearing bacteria. Those bacteria
# belong to the most-abundant species (i.e., species 1) if PReplMostAbun is
# TRUE, and to the least-abundant species (i.e., species nspecies) if
# PReplMostAbun is FALSE.
PReplMostAbun <- TRUE
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

## Parameter set to show that interactions do not influence the plots
intmeanset <- c(-1e-11, 5e-12)
selfintmeanset <- c(-1e-11, 0)
costset <- seq(from = 0, to = 0.1505, by = 0.0005)

## Parameter set to create detailed plot without showing effect of interactions
intmeanset <- 0
selfintmeanset <- -5e-12
costset <- seq(from = 0, to = 0.15025, by = 0.00025)


#### Running the simulations ####
set.seed(seed = 314, kind = "default", normal.kind = "default",
         sample.kind = "default")
starttime <- Sys.time()

# Create matrix to store data
# To do:
# - Create variable 'newgrowthratecode' as NA and let it be written to the data
#   such that the same number of columns are present for this and for PinNew.
nrowplotdata <- prod(lengths(list(nspeciesset, abunmodelset, intmeanset,
                                  selfintmeanset, costset, conjrateset,
                                  taxmattypeset),
                             use.names = FALSE))
print(paste(niter*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = 3*4 + 29)
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
      abundance <- brokenstick_fast(nspecies = nspecies, totalabun = totalabun,
                                    takelimit = TRUE)
      # Using a number instead of a name, to prevent type-conversion when
      # storing it in matrices.
      abunmodelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt_fast(nspecies = nspecies, totalabun = totalabun,
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
        data <- matrix(data = NA, nrow = nrowdata, ncol = 35)
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
            growthrate <- getgrowthrate(abundance = abundance, intmat = intmat)
            eqinfo <- geteqinfo(model = "gLV", abundance = abundance,
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
                  # initially plasmid-bearing species to the initially
                  # plasmid-free species. If 'PReplMostAbun' is TRUE the initially
                  # plasmid-bearing species is the most-abundant species (i.e.,
                  # species 1), such that the first row and column of the matrix
                  # taxmat have to be adjusted. If 'PReplMostAbun' is FALSE this
                  # is the least-abundant species (i.e., species nspecies), such
                  # that the last row and column of the matrix taxmat have to be
                  # adjusted. In both cases the diagonal should remain
                  # 'SameSpecies' because those reflect intraspecies
                  # relationships by definition.
                  
                  # Note that PReplMostAbun is used to determine the initially
                  # plasmid-bearing species is also the perturbed
                  # species.
                  if(PReplMostAbun == TRUE) {
                    taxmat[1, -1] <- taxmattype
                    taxmat[-1, 1] <- taxmattype
                  } else {
                    taxmat[nspecies, -nspecies] <- taxmattype
                    taxmat[-nspecies, nspecies] <- taxmattype
                  }
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
                                        abundance = c(abundance, rep(0, nspecies)),
                                        intmat = intmat, growthrate = growthrate,
                                        cost = cost, conjmat = conjmat)
                
                # Get equilibrium characteristics regarding ecological and
                # epidemiological stability of the plasmid-free equilibrium (see
                # Roberts and Heesterbeek 2021).
                eqinfoecol <- geteqinfo(model = "ecol",
                                        abundance = c(abundance, rep(0, nspecies)),
                                        intmat = intmat, growthrate = growthrate,
                                        cost = cost, conjmat = conjmat)
                eqinfoepi <- geteqinfo(model = "epi",
                                       abundance = c(abundance, rep(0, nspecies)),
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
                  cost, conjratecode, taxmatcode, iter, seq_len(nspecies),
                  abundance,
                  diag(intmat), c(growthrate), iterintmat,
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
                            "selfintdata", "growthrate",
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
                 summarydata))
        rowindexplotdata <- rowindexplotdatanew
      }
    }
  }
}
duration <- Sys.time() - starttime
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)

colnames(plotdata) <- c("niter", "nspecies", "abunmodelcode", "intmean",
                        "selfintmean",
                        colnames(summarydata))
summary(warnings())
rm(summarydata)


#### Saving settings and output to CSV and RDS files ####
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
if(PReplMostAbun == FALSE) {
  DateTimeStamp <- paste0(DateTimeStamp, "PReplLeastAbun")
}
names(conjrateset) <- paste0("conjrateset", seq_along(conjrateset))
settings <- c(list(niter = niter, niterintmat = niterintmat,
                   simulateinvasion = FALSE,
                   smallstate = smallstate, smallchange = smallchange,
                   tstep = NA,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   costset = costset, conjrateset, taxmattype = taxmattypeset,
                   costtype = costtype,
                   PFrom = if(PReplMostAbun) {"MostAbun"} else {"LeastAbun"},
                   PReplMostAbun = PReplMostAbun, duration = duration))
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

# Saving the data as R-object into an R data file takes much less space than
# saving it as csv. The R data files can be read into R using
# readRDS(file = file.path("OutputMS", "YYYY_MM_DD",
#                          "YYYY_MM_DD_MM_SS_plotdata.rds"))
saveRDS(object = plotdata,
        file = paste0(DateTimeStamp, "_plotdata.rds"))


#### Reading previously saved data from a CSV file ####
# # To read data from a CSV file, put the CSV file in the working directory
# # (see getwd()), uncomment this section and change the date-time-stamp in the
# # file name (the next line prints a list of files in the working directory
# # that contain 'multispecies' in their name).
# list.files(pattern = "multispecies", ignore.case = TRUE)
# # Note that the extension (.csv) should be included in the file name.
# filename <- file.path("OutputMS", "YYYY_MM_DD",
#                       "YYYY_MM_DD_MM_SSmultispecies.csv")
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
# DateTimeStamp <- "2025_03_28_15_58"
# plotdata <- readRDS(paste0(DateTimeStamp, "_plotdata.rds"))

# Create variables that do not yet exist if previous data was read into R
# instead of running the simulations.
if(!exists("conjratecode")) {
  conjratecode <- max(plotdata[, "conjratecode"])
}
if(!exists("abunmodel")) {
  abunmodel <- abunmodelset[length(abunmodelset)]
}


#### Labels and limits for plots ####
labspecies <- paste("Sp.", seq_len(maxnspecies))
names(labspecies) <- seq_len(maxnspecies)
labnspecies <- paste(nspeciesset, "sp.")
names(labnspecies) <- nspeciesset
labmodel <- c("Broken stick", "Dom. preemption")
names(labmodel) <- c(1, 2)
labcost <- paste0("Fitness cost\n", costset, "/h")
names(labcost) <- costset
labconjrate <- paste("Conjset", seq_along(conjrateset))
names(labconjrate) <- seq_along(conjrateset)
labtaxmat <- c("all\nconjugation\nrates equal",
               "lower\nintersp. conj.\nrates initP")
names(labtaxmat) <- seq_along(taxmattypeset)
# '.multi_line = FALSE' to collapse facet labels into a single label
mylabeller <- labeller(species = labspecies, nspecies = labnspecies,
                       abunmodelcode = labmodel,
                       cost = labcost, conjratecode = labconjrate,
                       taxmatcode = labtaxmat, .multi_line = FALSE,
                       .default = label_value)
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

## Show border of ecological stability with heatmap in CreatePlot_bifur()
# Values in a facet are either all stable, or all unstable, making it impossible
# to plot contours delimiting stable and unstable regions.
CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableecol",
           filltitle = "fracstableecol", filltype = "continuous", ratio = NULL,
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode + intmean", facety = "nspecies",
           filename = "ecostabxcostyconj")


## Show border of epidemiological stability with a contour plot in CreatePlot_bifur()
# 'save' is set to FALSE and ggsave() is used to ensure the added guides
# arguments are included in the saved plots.
CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
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

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
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

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
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

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
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

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
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
# Customise colour palette
val_name <- "taxmatcode"
p_nspecies <- as.factor(plotdata[, val_name, drop = TRUE])
# my_cols <- (scales::hue_pal()(nlevels(p_nspecies) + 1L))[
#   seq_len(nlevels(p_nspecies))]
my_cols <- c("#FFC20A", "#0C7BDC")
names(my_cols) <- levels(p_nspecies)

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
temp_plotdata$fracstableepi <- plotdata[, "conjrate"] * R1 + g22 * R2 + 
  sqrt(temp$p1) < 2 * plotdata[, "cost"]
plotdata_full <- rbind(cbind(plotdata, type = "numeric"),
                       cbind(temp_plotdata, type = "analytic"))
rm(row_ind)
rm(temp)
rm(temp_plotdata)

conjrate_inter_reduced <- getconjmat(
  nspecies = 2, conjrate = rep(conjrate_base, 2L),
  taxmat = matrix(c("SameSpecies", "OtherClass", "OtherClass", "SameSpecies"),
                  nrow = 2, byrow = FALSE))[1, ]

plotdata_full$type <- factor(plotdata_full$type, levels = c("numeric", "analytic"))

p_comp_epistab <- CreatePlot_bifur(
  dataplot = plotdata_full,
  xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
  contour_var = "fracstableepi", contour_col = "taxmatcode", contour_lty = "type",
  limx = range(c(0, costset)), limy = range(log10(seqconjrate)), ratio = NULL,
  labx = "Fitness cost of bearing a plasmid",
  laby = paste0("Log10(intraspecies conjugation rate",
                " of\nthe initially plasmid-bearing",
                " species)"),
  linezero = FALSE, facetx = ".", facety = "nspecies", rotate_x_labels = FALSE,
  save = FALSE) +
  theme(legend.box = "horizontal",
        legend.direction = "vertical",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt"),
        legend.title = element_blank()) +
  scale_colour_manual(values = my_cols,
                      labels = c("All conjugation rates equal",
                                 paste0("Lower interspecies conjugation rates",
                                        " to and\nfrom the initially",
                                        " plasmid-bearing species"))) +
  labs(caption = NULL)
p_comp_epistab
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabynspecies_comparison.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

plotdata_2sp <- plotdata_full[plotdata_full$nspecies == 2, , drop = FALSE]
costmark_2sp <- rep((conjrate_base * abundance[2] +
                       sqrt((-conjrate_base * abundance[2])^2 +
                              4 * (conjrate_inter_reduced)^2 * abundance[1] * abundance[2]))/2,
                    each = nrow(plotdata_2sp) / length(conjrate_inter_reduced))

p_comp_epistab +
  geom_vline(data = plotdata_2sp, aes(xintercept = costmark_2sp),
             show.legend = FALSE, linetype = 2) +
  guides(linetype = "none")
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabynspecies_comparison_nolegend_withlines.png"),
         width = 2150 / 1.27, height = 2150, units = "px", dpi = 300)
  ggsave(paste0(DateTimeStamp, "FigS06.png"),
         width = 2150 / 1.27, height = 2150, units = "px", dpi = 300)
}
rm(costmark_2sp)
rm(plotdata_2sp)
rm(p_comp_epistab)

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "intmean",
           contour_lty = "selfintmean",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate",
                         " of\nthe initially plasmid-bearing",
                         " species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspeciescolltyinter.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

# Note: assuming sets are chosen in such a way that border of invasion is shown
# in the plot
nspecies_numeric <- as.numeric(as.character(plotdata[, "nspecies"]))
plotdata_manysp <- plotdata[nspecies_numeric > max(nspecies_numeric) - 0.25,
                            , drop = FALSE]
min_conj_rate <- min(plotdata_manysp[, "conjrate"], na.rm = TRUE)
col_split <- c("conjratecode", "taxmatcode")
bool_ind_asymp_stableepi <- plotdata_manysp[, "conjrate"] == min_conj_rate &
  plotdata_manysp[, "fracstableepi"] == 1
data_asymp_stableepi <- split(x = plotdata_manysp[bool_ind_asymp_stableepi, "cost"],
                              f = plotdata_manysp[bool_ind_asymp_stableepi, col_split],
                              drop = TRUE, sep = "_")

bool_ind_asymp_unstableepi <- plotdata_manysp[, "conjrate"] == min_conj_rate &
  plotdata_manysp[, "fracstableepi"] == 0
data_asymp_unstableepi <- split(x = plotdata_manysp[bool_ind_asymp_unstableepi, "cost"],
                                f = plotdata_manysp[bool_ind_asymp_unstableepi, col_split],
                                drop = TRUE, sep = "_")

asymp_cost <- colMeans(rbind(unlist(lapply(X = data_asymp_stableepi,
                                           FUN = min, na.rm = TRUE)),
                             unlist(lapply(X = data_asymp_unstableepi,
                                           FUN = max, na.rm = TRUE))))
offset_x <- (min(asymp_cost) - min(costset)) / 2.4 # / 4
offset_y <- (log10(max(plotdata[, "conjrate"])) -
               log10(max(plotdata[plotdata[, "fracstableepi"] > 0.9, "conjrate"]))) / 4
p_comp_epistab_v2 <- CreatePlot_bifur(dataplot = plotdata_manysp,
                                xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
                                contour_var = "fracstableepi", contour_col = val_name,
                                contour_lty = val_name,
                                limx = range(c(0, costset)),
                                limy = range(log10(seqconjrate)),
                                ratio = NULL,
                                labx = "Fitness cost of bearing a plasmid",
                                laby = paste0("Log10(intraspecies conjugation rate",
                                              " of\nthe initially plasmid-bearing",
                                              " species)"),
                                linezero = FALSE, facetx = ".", facety = ".",
                                base_size = 17, rotate_x_labels = FALSE,
                                save = FALSE) +
  theme(legend.box.margin = margin(-5, 0, 0, -75),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = "left",
        legend.margin = margin(-5, 0, -5, 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = margin(10, 12, 10, 5)) +
  guides(linetype = "none") +
  scale_colour_manual(values = c('1' = "#FFC20A", '2' = "#0C7BDC"),
                      labels = c("All conjugation rates equal",
                                 paste0("Lower interspecies conjugation rates",
                                        " to and\nfrom the initially",
                                        " plasmid-bearing species"))) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2) +
  geom_vline(xintercept = asymp_cost, show.legend = FALSE, linetype = 2) +
  labs(caption = NULL) +
  annotate("label",
           label = paste0("Plasmid invasion possible even if\ninterspecies",
                          " conjugation rates to\nand from initially",
                          " plasmid-bearing\nspecies are decreased"),
           x = min(costset) + offset_x,
           y = log10(max(plotdata[, "conjrate"])) - offset_y,
           hjust = "inward", vjust = "inward", size = 4) +
  annotate("label",
           label = paste0("Plasmid invasion not possible even if",
                          "\ninterspecies conjugation rates to\nand from",
                          " initially plasmid-bearing\nspecies are not decreased"),
           x = max(costset) - offset_x,
           y = log10(min(plotdata[, "conjrate"])) + offset_y,
           hjust = "inward", vjust = "inward", size = 4)
p_comp_epistab_v2
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistab_v2a.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}
rm(plotdata_manysp)

probs_y <- 0.73
p_comp_epistab_v2 +
  annotate("text",
           label = c("1a", "2a", "2b", "3a", "3b", "3c"),
           x = c(mean(c(min(costset), min(asymp_cost))),
                 mean(asymp_cost),
                 mean(asymp_cost),
                 mean(c(max(costset), max(asymp_cost))),
                 mean(c(max(costset), max(asymp_cost))),
                 mean(c(max(costset), max(asymp_cost)))),
           y = c(mean(log10(range(plotdata[, "conjrate"], na.rm = TRUE))),
                 quantile(log10(range(plotdata[, "conjrate"], na.rm = TRUE)),
                          probs = probs_y),
                 quantile(log10(range(plotdata[, "conjrate"], na.rm = TRUE)),
                          probs = 1 - probs_y),
                 quantile(log10(range(plotdata[, "conjrate"], na.rm = TRUE)),
                          probs = mean(c(probs_y, 1))),
                 quantile(log10(range(plotdata[, "conjrate"], na.rm = TRUE)),
                          probs = probs_y),
                 quantile(log10(range(plotdata[, "conjrate"], na.rm = TRUE)),
                          probs = 1 - probs_y)),
           hjust = "inward", vjust = "inward", size = 4)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistab_v2b.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
  ggsave(paste0(DateTimeStamp, "Fig02.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

# Need to set filltype to continuous to prevent error on missing filllabels
CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableepi",
           filltitle = "fracstableepi", contour_var = NULL, contour_col = NULL,
           contour_lty = NULL, filltype = "continuous", ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate",
                         " of\nthe initially plasmid-bearing",
                         " species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
           rotate_legend = TRUE,
           filename = "epistabheatmap")

# To do:
# - Could also use this to identify the coordinates of the intersections and
#   draw the relevant line segments and place labels in the correct positions
#   instead of filling the parts with colour.
data_taxmatcode_1 <- plotdata[plotdata$taxmatcode == 1, ]
data_taxmatcode_1$data_comp_epi <- NA
data_taxmatcode_1$data_comp_epi_char <- NA_character_

ind_both_zero <- which(data_taxmatcode_1$fracstableepi == 0 &
                         plotdata[plotdata$taxmatcode == 2, "fracstableepi"] == 0)
data_taxmatcode_1$data_comp_epi[ind_both_zero] <- 0
data_taxmatcode_1$data_comp_epi_char[ind_both_zero] <- "inv"

ind_intermediate <- which(plotdata[plotdata$taxmatcode == 1, "fracstableepi"] !=
                            plotdata[plotdata$taxmatcode == 2, "fracstableepi"])
data_taxmatcode_1$data_comp_epi[ind_intermediate] <- 0.5
data_taxmatcode_1$data_comp_epi_char[ind_intermediate] <- "someinv"

ind_both_one <- which(data_taxmatcode_1$fracstableepi == 1 &
                        plotdata[plotdata$taxmatcode == 2, "fracstableepi"] == 1)
data_taxmatcode_1$data_comp_epi[ind_both_one] <- 1
data_taxmatcode_1$data_comp_epi_char[ind_both_one] <- "noinv"

data_taxmatcode_1$data_comp_epi_fac <- factor(data_taxmatcode_1$data_comp_epi,
                                              levels = c(0, 0.5, 1))
levels(data_taxmatcode_1$data_comp_epi_fac) <- c("allinv", "intinv", "noinv")

data_taxmatcode_1$data_comp_cost <- NA
data_taxmatcode_1$data_comp_cost[which(data_taxmatcode_1$cost < min(asymp_cost))] <- "lowcost"
data_taxmatcode_1$data_comp_cost[which(data_taxmatcode_1$cost >= min(asymp_cost) &
                                         data_taxmatcode_1$cost <= max(asymp_cost))] <- "intcost"
data_taxmatcode_1$data_comp_cost[which(data_taxmatcode_1$cost > max(asymp_cost))] <- "highcost"

data_taxmatcode_1$data_comp_epi_cost <- NA
data_taxmatcode_1$data_comp_epi_cost <- interaction(data_taxmatcode_1$data_comp_epi_char,
                                                    data_taxmatcode_1$data_comp_cost,
                                                    drop = TRUE)

CreatePlot_bifur(dataplot = data_taxmatcode_1,
           xvar = "cost", yvar = "log10(conjrate)",
           fillvar = "as.numeric(data_comp_epi_cost)",
           filltitle = "", contour_var = NULL, contour_col = NULL,
           contour_lty = NULL, filltype = "continuous", ratio = NULL,
           title = paste0("Effect of reducing interspecies\nconjugation rate on",
                          " plasmid invasion"),
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate",
                         " of\nthe initially plasmid-bearing",
                         " species)"),
           linezero = FALSE, facetx = ".", facety = "nspecies",
           rotate_x_labels = FALSE, save = FALSE) +
  theme(legend.position = "none") + labs(caption = NULL)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "compareinvasion.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

# Need to be as many labels as there are levels in fillvar
levels(data_taxmatcode_1$data_comp_epi_cost)
filllabels <- letters[seq_len(nlevels(data_taxmatcode_1$data_comp_epi_cost))]
names(filllabels) <- filllabels

levels(data_taxmatcode_1$data_comp_epi_cost)
filllabels <- levels(data_taxmatcode_1$data_comp_epi_cost)
names(filllabels) <- filllabels

CreatePlot_bifur(dataplot = data_taxmatcode_1,
           xvar = "cost", yvar = "log10(conjrate)",
           fillvar = "data_comp_epi_cost",
           filltitle = "Region", contour_var = NULL, contour_col = NULL,
           contour_lty = NULL, filltype = "discrete", limits = filllabels,
           ratio = NULL,
           title = paste0("Effect of reducing interspecies\nconjugation rate on",
                          " plasmid invasion"),
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate",
                         " of\nthe initially plasmid-bearing",
                         " species)"),
           linezero = FALSE, facetx = ".", facety = "nspecies",
           rotate_x_labels = FALSE, save = FALSE) +
  labs(caption = NULL) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE))
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "compareinvasion_withlegend.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}
