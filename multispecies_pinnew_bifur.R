################################################################################
## Modelling the effects of ecological interactions and distinct conjugation  ##
## rates on the invasion of a conjugative plasmid in bacterial communities    ##
################################################################################


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
source("./multispecies_funcs.R")


#### Settings and defining parameter space ####

##### Basis parameter set #####
totalabun <- 1e11
abunmodelset <- c("dompreempt")
nspeciesset <- 1 + c(2, 4, 8, 16) # Because species do matter for threshold of costs
maxnspecies <- max(nspeciesset)
newgrowthratecode <- 2
costmark <- NULL # Plot dotted vertical lines at indicated values if not NULL
# The conjugation rate given here is the 'overall' conjugation rate. For all
# species, the intraspecies conjugation rate will be equal to this overall
# conjugation rate. The interspecies conjugation rate will be equal to this
# overall conjugation rate if all populations belong to the same species.
conjrate_base <- 1e-12
seqconjrate <- 10^seq(from = -13.0, to = -11.5, by = 0.05)
conjrateset <- NULL
for(conjrate in seqconjrate) {
  conjrateset <- c(conjrateset, list(rep(conjrate, maxnspecies)))
}
# If taxmattype is "SameSpecies", the conjugation rate is the same for all
# populations, and equal to 'conjrateset' given above. If taxmattype is
# "OtherClass", the interspecies conjugation rate to and from the initially
# plasmid-bearing population (the newly added species 1) is reduced a 1000-fold.
taxmattypeset <- c("SameSpecies", "OtherClass")
smallstate <- NA
finalsmallstate <- NA
smallchange <- NA
saveplots <- TRUE
niterintmat <- 1
nsimul <- 1
intmeanset <- c(-1e-11, 0, 5e-12)
selfintmeanset <- c(-1e-11, -5e-12, 0)
costset <- seq(from = 0, to = 1, by = 0.0025)
# Repeat 8 colours to plot 16 species.
mycol <- rep(c("black", "blue", "red", "darkgreen", "darkgrey", "brown", "purple",
               "darkorange", "green1", "yellow", "hotpink"), 2)


#### Running the simulations ####
set.seed(seed = 314, kind = "default", normal.kind = "default", sample.kind = "default")
starttime <- Sys.time()

# Create matrix to store data
# To do:
# - Create variable 'newgrowthratecode' as NA and let it be written to the data
#   such that the same number of columns are present for this and for PinNew.
nrowplotdata <- prod(lengths(list(nspeciesset, abunmodelset, intmeanset,
                                  selfintmeanset, costset, conjrateset, taxmattypeset),
                             use.names = FALSE))
print(paste(nsimul*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = 3*4 + 30)
nrowdatatotal <- prod(lengths(list(abunmodelset,intmeanset, selfintmeanset,
                                   costset, conjrateset, taxmattypeset),
                              use.names = FALSE))*nsimul*sum(nspeciesset)
indexdatatotal <- 1

# Run simulations
rowindexplotdata <- 1
rowindexdata <- 1

for(nspecies in nspeciesset) {
  for(abunmodel in abunmodelset) {
    if(abunmodel == "brokenstick") {
      abundance <- brokenstick_fast(nspecies = nspecies - 1, totalabun = totalabun,
                                    takelimit = TRUE)
      # Using a number instead of a name, to prevent type-conversion when
      # storing it in matrices.
      abunmodelcode <- 1
    }
    if(abunmodel == "dompreempt") {
      abundance <- dompreempt_fast(nspecies = nspecies - 1, totalabun = totalabun,
                                   takelimit = TRUE)
      abunmodelcode <- 2
    }
    
    for(intmean in intmeanset) {
      
      for(selfintmean in selfintmeanset) {
        print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
                     ", intmean = ", intmean, ", selfintmean = ", selfintmean,
                     ": started at ", Sys.time()), quote = FALSE)
        nrowdata <- nsimul * nspecies * length(costset) * length(conjrateset) *
          length(taxmattypeset)
        data <- matrix(data = NA, nrow = nrowdata, ncol = 36)
        indexdata <- 1
        relabunRsp <- rep(NA, maxnspecies)
        relabunRconjsp <- rep(NA, maxnspecies)
        relabunPconjsp <- rep(NA, maxnspecies)
        
        for(iter in seq_len(nsimul)) {
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
            growthrateeq <- getgrowthrate(abundance = abundance,
                                          intmat = intmat[-1, -1])
            # The growth rate of new species is the mean growth rate of the
            # plasmid-free species decreased with 2 standard deviations,
            # unchanged, or increased with 2 standard deviations when
            # 'newgrowthratecode' is 1, 2, or 3, respectively.
            growthrate <- c(switch(newgrowthratecode,
                                   mean(growthrateeq) - 2*sd(growthrateeq),
                                   mean(growthrateeq),
                                   mean(growthrateeq) + 2*sd(growthrateeq)),
                            growthrateeq)
            eqinfo <- geteqinfo(model = "gLV", abundance = c(0, abundance),
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
                  # initially plasmid-bearing species to the initially plasmid-
                  # free species. In this model the initially plasmid-bearing
                  # species is always the newly-introduced species 1, such that
                  # the first row and column of the matrix taxmat have to be
                  # adjusted. The diagonal should remain 'SameSpecies' because
                  # those reflect intraspecies relationships by definition.
                  taxmat[1, -1] <- taxmattype
                  taxmat[-1, 1] <- taxmattype
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
                                        abundance = c(0, abundance, rep(0, nspecies)),
                                        intmat = intmat, growthrate = growthrate,
                                        cost = cost, conjmat = conjmat)
                
                # Get equilibrium characteristics regarding ecological and
                # epidemiological stability of the plasmid-free equilibrium (see
                # Roberts and Heesterbeek 2021). The abundance of 0 is the
                # initially plasmid-bearing population.
                eqinfoecol <- geteqinfo(model = "ecol",
                                        abundance = c(0, abundance, rep(0, nspecies)),
                                        intmat = intmat, growthrate = growthrate,
                                        cost = cost, conjmat = conjmat)
                eqinfoepi <- geteqinfo(model = "epi",
                                       abundance = c(0, abundance, rep(0, nspecies)),
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
                  nsimul, nspecies, abunmodelcode, intmean, selfintmean,
                  cost, conjratecode, taxmatcode, iter, seq_len(nspecies), c(0, abundance),
                  diag(intmat), c(growthrate), newgrowthratecode, iterintmat,
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
        
        colnames(data) <- c("nsimul", "nspecies", "abunmodelcode",
                            "intmean", "selfintmean", "cost", "conjratecode",
                            "taxmatcode", "iter", "species", "abundance",
                            "selfintdata", "growthrate", "newgrowthratecode",
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
          tibble(nsimul, nspecies, abunmodelcode, intmean, selfintmean,
                 newgrowthratecode, summarydata))
        rowindexplotdata <- rowindexplotdatanew
      }
    }
  }
}
duration <- Sys.time() - starttime
print(paste0("Finished simulations: ", Sys.time()), quote = FALSE)

colnames(plotdata) <- c("nsimul", "nspecies", "abunmodelcode", "intmean",
                        "selfintmean", "newgrowthratecode",
                        colnames(summarydata))
summary(warnings())
rm(summarydata)


#### Saving settings and output to CSV and RDS files ####
DateTimeStamp <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M"), "PInNewSp")
names(conjrateset) <- paste0("conjrateset", seq_along(conjrateset))
settings <- c(list(nsimul = nsimul, niterintmat = niterintmat,
                   simulateinvasion = FALSE,
                   smallstate = smallstate, smallchange = smallchange,
                   tstep = NA,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   newgrowthratecode = newgrowthratecode,
                   costset = costset, conjrateset, taxmattype = taxmattypeset,
                   PFrom = "PInNewSp", PReplMostAbun = TRUE, duration = duration))
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
#                       "YYYY_MM_DD_MM_SSPInNewSpmultispecies.csv")
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


#### Labels and limits for plots ####
labspecies <- paste("Sp.", seq_len(maxnspecies))
names(labspecies) <- seq_len(maxnspecies)
labnspecies <- paste(nspeciesset, "sp.")
names(labnspecies) <- nspeciesset
abunmodel_lab <- c("Broken stick", "Dom. preemption")
names(abunmodel_lab) <- c(1, 2)
labcost <- paste0("Fitness cost\n", costset, "/h")
names(labcost) <- costset
labconjrate <- paste("Conjset", seq_along(conjrateset))
names(labconjrate) <- seq_along(conjrateset)
labtaxmat <- c("all conj. rates\nequal",
               "lower intersp.\nconj. rates initP\n")
names(labtaxmat) <- seq_along(taxmattypeset)
# '.multi_line = FALSE' to collapse facet labels into a single label
mylabeller <- labeller(species = labspecies, nspecies = labnspecies,
                       abunmodelcode = abunmodel_lab,
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
           rotate_legend = TRUE,
           filename = "ecostabxcostyconj", width = 1.75 * 1850, height = 2675)


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
         width = 1.75 * 1850, height = 2675, units = "px", dpi = 300)
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
  theme(legend.box = "vertical",
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
temp_plotdata$fracstableepi <- plotdata[, "conjrate"] * R1 + g22 * R2 + sqrt(temp$p1) <
  2 * plotdata[, "cost"]
rm(temp)
rm(row_ind)

temp_plotdata <- rbind(cbind(plotdata, type = "numeric"),
                       cbind(temp_plotdata, type = "analytic"))
CreatePlot_bifur(dataplot = temp_plotdata,
           xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "as.factor(type)",
           contour_lty = "taxmatcode",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Comparing epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = ".", facety = "nspecies",
           rotate_x_labels = FALSE, save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabynspecies_comparison.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}
rm(temp_plotdata)

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "intmean",
           contour_lty = "selfintmean",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode", facety = "nspecies",
           save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt"),
        legend.position = "bottom",
        panel.border = element_blank(),
        panel.spacing.x = unit(3, "pt"),
        panel.spacing.y = unit(6, "pt"),
        strip.background = element_rect(color = NA)) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  labs(caption = NULL) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistabxtaxmatynspeciescolltyinter.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
  ggsave(paste0(DateTimeStamp, "Fig20.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

# Note: assuming sets are chosen such that border of invasion is shown in the
# plot
CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           contour_var = "fracstableepi", contour_col = "nspecies",
           contour_lty = "taxmatcode",
           limx = range(c(0, costset)), limy = range(log10(seqconjrate)),
           ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = ".", facety = ".",
           rotate_x_labels = FALSE, save = FALSE) +
  theme(legend.box = "horizontal",
        legend.margin = margin(c(-5, 0, -5, 0), unit = "pt")) +
  guides(col = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  geom_vline(xintercept = costmark, show.legend = FALSE, linetype = 2) +
  annotate("text", label = c("Invasion of plasmid-\nbearing bacteria\n[WARNING: should include (self)intmean!",
                             "No invasion of plasmid-\nbearing bacteria"),
           x = quantile(costset, c(0.1, 0.9)),
           y = quantile(log10(seqconjrate), c(0.9, 0.1)), 
           hjust = "inward", vjust = "inward", size = 4)
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "epistab.png"),
         width = 2150, height = 2150, units = "px", dpi = 300)
}

# Need to set filltype to continuous to prevent error on missing filllabels
CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = NULL,
           filltitle = NULL, contour_var = "fracstableepi", contour_col = NULL,
           contour_lty = NULL, filltype = "continuous", ratio = NULL,
           title = "Epidemiological (in)stability",
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "intmean + selfintmean",
           facety = "nspecies + taxmatcode",
           rotate_x_labels = TRUE, save = TRUE)

CreatePlot_bifur(xvar = "cost", yvar = "log10(conjrate)", fillvar = "fracstableepi",
           filltitle = "fracstableepi", contour_var = NULL, contour_col = NULL,
           contour_lty = NULL, filltype = "continuous", ratio = NULL,
           title = NULL,
           labx = "Fitness cost of bearing a plasmid",
           laby = paste0("Log10(intraspecies conjugation rate of\nthe",
                         " initially plasmid-bearing species)"),
           linezero = FALSE, facetx = "taxmatcode + intmean + selfintmean",
           facety = "nspecies", rotate_x_labels = TRUE, width = 9*1650/2,
           height = 2675, filename = "epistabheatmap_morefacets_v3")
# NB. The conjugation rate of the initially plasmid-bearing species does not
# have an effect, and the fitness costs of bearing a plasmid has only a single
# relevant value for each combination of nspecies * taxmatcode * intmean *
# selfintmean. So it would be much more interesting to have the save layout as
# other figures in the main text (e.g., Fig. 3) but with the heatmap indicating
# the minimum costs at which the plasmid can invade! So change the dataset to a
# larger range of interactions, a single conjugation rate (i.e., the default
# 1e-12), and a range of costs (note that costs are not relative, so using a
# range from 0 to 1 for costs is confusing!
