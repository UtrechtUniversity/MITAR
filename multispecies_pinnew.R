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


#### Optionally to do ####
# The matrix 'data' is converted to a tibble to efficiently get summary
# statistics. The tibble with summary statistics is then converted to a matrix
# to fill in part of the matrix plotdata. Using a tibble for 'data' from the
# start prevents some of this type-conversion, might make naming columns on
# assignment clearer, and enable use of the more efficient dplyr::bind_cols()
# instead of base::cbind(). However, originally using a data.frame instead of
# matrix to store data made progress really slow, so I should check how timing
# is affected when 'data', 'plotdata' or both are tibbles from the start.


#### Loading required packages ####
library(deSolve)   # checkequilibrium() and perturbequilibrium_pinnew() call ode()
                   # if showplot and simulateinvasion are 'TRUE', respectively
library(dplyr)     # across(), full_join(), group_by(), near(), summarise()
library(ggplot2)   # to display data and results
library(scales)    # CreatePlot() calls format_scaled() or format_sci() which
                   # call scales::label_number() or scales::label_scientific(),
                   # respectively, to format values on the axes.
library(rootSolve) # geteqinfo() calls jacobian.full()
library(TruncatedNormal) # getintmat calls rtnorm()
# On the pipe operator (%>%), see ?'%>%' and Ch. 18 'Pipes' in Wickham 2017.


#### Load required functions ####
source("./multispecies_funcs.R")


#### Settings and defining parameter space ####
##### Notes #####
# - See the annotation of the functions that use the arguments for more detailed
#   information.
# - To simulate that each species belongs to a different class, an additional
#   taxmattype has to be added. Then only the interspecies conjugation rates
#   should be reduced, as the intraspecies conjugation rates of all species are
#   still those given in conjrateset. This unchanged intraspecies conjugation
#   rate makes it different from reducing conjrateset 1000-fold and using
#   'taxmatsame' as taxmattype.

##### Basis parameter set #####
totalabun <- 1e11
abunmodel_options <- c("brokenstick", "dompreempt", "equalabun")
abunmodel_lab <- c("Broken stick", "Dom. preemption", "Equal abundance")
abunmodelset <- abunmodel_options[1] # Can also use multiple values
# The growth rate of new species is the mean growth rate of the plasmid-free
# species decreased with 2 standard deviations, unchanged, or increased with 2
# standard deviations when newgrowthratecode == 1, 2, or 3, respectively.
newgrowthratecode <- 2 # A SINGLE value (1, 2 or 3) should be provided.
costset <- c(0.05, 0.09)
# The conjugation rate given here is the 'overall' conjugation rate. For all
# species, the intraspecies conjugation rate will be equal to this overall
# conjugation rate. See 'taxmattypeset' on the interspecies conjugation rates.
conjrate_base <- 1e-12
# If taxmattype is "SameSpecies", the conjugation rate is the same for all
# populations, and equal to 'conjrateset' given above. If taxmattype is
# "OtherClass", the interspecies conjugation rate to and from the initially
# plasmid-bearing population (the newly added species 1) is reduced a 1000-fold.
taxmattypeset <- c("SameSpecies", "OtherClass")
# Some plasmid-free bacteria are added to simulate perturbation by plasmid-free
# bacteria, and some plasmid-free bacteria are replaced with plasmid-bearing
# bacteria to simulate perturbation by plasmid-bearing bacteria. Those bacteria
# belong to the most-abundant species (i.e., species 2, because species 1 is the
# newly-introduced species) if PReplMostAbun is TRUE, and to the least-abundant
# species (i.e., species nspecies) if PReplMostAbun is FALSE.
PReplMostAbun <- TRUE
simulateinvasion <- TRUE # Should simulations over time be performed?
# Variables that become smaller than smallstate during the integration are set to 0.
smallstate <- 1e-3
finalsmallstate <- 1
# If the sum of absolute rates of change is equal to smallchange, equilibrium is
# assumed to be reached and the integration is terminated.
smallchange <- 1e-2
niterintmat <- 1
saveplots <- TRUE
nspeciesset <- 1 + c(2, 4, 8, 16)
intmeanset <- seq(from = -1e-11, to = 5e-12, by = 5e-13)
selfintmeanset <- seq(from = -1e-11, to = 0, by = 5e-13)
nsimul <- 50
saveplotconjovertime <- FALSE
# Repeat 8 colours to plot 16 species.
mycol <- rep(c("black", "blue", "red", "darkgrey", "darkorange", "green1",
               "yellow", "hotpink"), 2)

##### Parameter set to show plots over time (FigS08) #####
# nspeciesset <- 1 + 4
# intmeanset <- 0
# selfintmeanset <- -0.5e-11
# nsimul <- 1
# saveplotconjovertime <- TRUE
# mycol <- rep(c("#000000", "#E41A1C", "#377EB8", "#FF7F00", "#999999"), 2)


#### Running the simulations ####
maxnspecies <- max(nspeciesset)
names(abunmodel_lab) <- seq_along(abunmodel_lab)
conjrateset <- list(rep(conjrate_base, maxnspecies))
set.seed(seed = 314, kind = "default", normal.kind = "default",
         sample.kind = "default")
starttime <- Sys.time()

# Create matrix to store data
nrowplotdata <- prod(lengths(list(nspeciesset, abunmodelset, intmeanset,
                                  selfintmeanset, costset, conjrateset,
                                  taxmattypeset),
                             use.names = FALSE))
# When changing the number of columns, also update the line
# 'data <- matrix(data = NA, nrow = nrowdata, ncol = 60 + 4*maxnspecies)' below.
ncolplotdata <- if(simulateinvasion == TRUE) {
  18*4 + 4*4*maxnspecies + 39
} else {
  3*4 + 30
}
print(paste(nsimul*nrowplotdata, "simulations to run."), quote = FALSE)
plotdata <- matrix(data = NA, nrow = nrowplotdata, ncol = ncolplotdata)
nrowdatatotal <- prod(lengths(list(abunmodelset,intmeanset, selfintmeanset,
                                   costset, conjrateset, taxmattypeset),
                              use.names = FALSE))*nsimul*sum(nspeciesset)
datatotal <- matrix(data = NA, nrow = nrowdatatotal, ncol = 60 + 4*maxnspecies)
indexdatatotal <- 1

# Run simulations
rowindexplotdata <- 1
rowindexdata <- 1

for(nspecies in nspeciesset) {
  # Note: pertpop and pertpopconj indicate which population is perturbed by
  # adding bacteria, which in this model is always the newly-introduced species
  # 1. In case of pertpopconj the added bacteria are plasmid-bearing and replace
  # plasmid-free bacteria. In case of pertpop, the bacteria are plasmid-free and
  # are added to the existing bacteria.
  pertpop <- "R1"
  pertpopconj <- "P1"
  
  for(abunmodel in abunmodelset) {
    # Use nspecies - 1 because the newly-introduced species (species 1) should
    # not yet be taken into account.
    abundance <- switch(
      abunmodel,
      brokenstick = brokenstick_fast(nspecies = nspecies - 1,
                                     totalabun = totalabun, takelimit = TRUE),
      dompreempt = dompreempt_fast(nspecies = nspecies - 1,
                                   totalabun = totalabun, takelimit = TRUE),
      equalabun = equalabun(nspecies = nspecies - 1, totalabun = totalabun),
      stop("'abunmodel' should be one of '",
           paste(abunmodel_options, collapse = "', '"),
           "', not '", abunmodel, "'."))
    
    # Using a number instead of a name, to prevent type-conversion when
    # storing it in matrices.
    abunmodelcode <- match(x = abunmodel, table = abunmodel_options)
    
    for(intmean in intmeanset) {
      
      for(selfintmean in selfintmeanset) {
        print(paste0("nspecies = ", nspecies, ", abundance model = ", abunmodel,
                     ", intmean = ", intmean, ", selfintmean = ", selfintmean,
                     ": started at ", Sys.time()), quote = FALSE)
        nrowdata <- nsimul * nspecies * length(costset) * length(conjrateset) *
          length(taxmattypeset)
        data <- matrix(data = NA, nrow = nrowdata, ncol = 60 + 4*maxnspecies)
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
            # The growth rate of the new species is the mean growth rate of the
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
          
          # Simulate invasion of plasmid-free bacteria into the plasmid-free
          # equilibrium in the model without conjugation (i.e., test internal
          # stability).
          if(iter == 1L) {
            warning("'tmax' is put at 100 for pertubation with plasmid-free",
                    " bacteria, such that\nequilibrium will be rarely reached.",
                    " This is done to save time and focus on\ninvasion of",
                    " plasmid-bearing bacteria.")
          }
          if(simulateinvasion == TRUE) {
            if(stableeq == FALSE) {
              abunfinal <- perturbequilibrium_pinnew(abundance = c(0, abundance),
                                              intmat = intmat,
                                              growthrate = growthrate,
                                              cost = NULL, conjmat = NULL,
                                              model = "gLV", pertpop = pertpop,
                                              # HARDCODED tmax = 100 to focus on
                                              # invasion of plasmid-bearing
                                              # bacteria, see the warning issued
                                              # for iter == 1L above.
                                              pertmagn = 1, tmax = 100, tstep = 10,
                                              showplot = FALSE, verbose = FALSE,
                                              silentinfgrowth = TRUE,
                                              silenteqnotreached = TRUE)
            } else {
              # No need for simulations if equilibrium is stable. Use nspecies
              # - 1 because the newly-introduced species (species 1) will go
              # extinct and thus should not be counted.
              abunfinal <- list(R = c(0, abundance), Rtotal = sum(abundance),
                                npopR = nspecies - 1,
                                P = NULL, Ptotal = NULL, npopP = NULL,
                                pertpopconjsurvived = NULL,
                                timepertpopconjextinct = NULL,
                                timefinal = 1, tmaxshort = 0,
                                eqreached = 1, infgrowth = 0)
            }
          } else {
            # No simulations over time performed, so set values to NA
            abunfinal <- list(R = rep(NA, nspecies), Rtotal = NA,
                              npopR = NA,
                              P = NULL, Ptotal = NULL, npopP = NULL,
                              pertpopconjsurvived = NULL,
                              timepertpopconjextinct = NULL,
                              timefinal = NA, tmaxshort = NA,
                              eqreached = NA, infgrowth = NA)
          }
          
          # Model without conjugation, so Rtotal is total abundance because P
          # does not exist
          relabunRsp[seq_len(nspecies)] <- abunfinal$R / abunfinal$Rtotal
          
          for(cost in costset) {
            conjratecode <- 0
            
            # Use conjrate <- conjrateset[[1]] to manually select the first
            # element
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
                  # those values reflect intraspecies relationships by definition.
                  taxmat[1, -1] <- taxmattype
                  taxmat[-1, 1] <- taxmattype
                }
                conjmat <- getconjmat_fast(nspecies = nspecies,
                                      conjrate = conjrate, taxmat = taxmat)
                
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
                
                # Append the signs of the real parts of the largest eigenvalues
                # to indicate change in stability without and with the plasmid
                # (-1 = stable, 0 = neutral, 1 = unstable). The + 3 ensure
                # numbers are positive, to prevent paste0(-1, -1) leading to NA
                # because of as.integer(-1-1))
                compstability <- as.integer(paste0(eqinfo["eigvalReSign"] + 3,
                                                   eqinfoconj["eigvalReSign"] + 3))
                
                # To simulate invasion of plasmid-bearing bacteria into the
                # plasmid-free equilibrium, the abundances of the plasmid-free
                # populations have to be appended to the abundances of the
                # plasmid-bearing populations.
                if(simulateinvasion == TRUE) {
                  if(eqinfoconj["eigvalRe"] >= 0 || saveplotconjovertime) {
                    abunfinalconj <- perturbequilibrium_pinnew(
                      abundance = c(0, abundance, rep(0, nspecies)),
                      intmat = intmat, growthrate = growthrate,
                      cost = cost, conjmat = conjmat,
                      model = "gLVConj", pertpop = pertpopconj,
                      pertmagn = 1, tmax = 5e3, tstep = 10,
                      showplot = saveplotconjovertime, ylim = c(1e5, 1e11),
                      addline = FALSE, verbose = FALSE,
                      silentinfgrowth = TRUE,
                      silenteqnotreached = TRUE)
                  } else {
                    # No need for simulations if equilibrium is stable
                    # NOTE: timepertpopconjextinct = 0 is only true if pertmagn
                    # is equal to finalsmallstate. Similarly, timefinal = 1 is
                    # not strictly true, as it will take some time to reach the
                    # new equilibrium.
                    abunfinalconj <- list(R = c(0, abundance),
                                          Rtotal = sum(abundance),
                                          npopR = nspecies - 1,
                                          P = rep(0, nspecies), Ptotal = 0,
                                          npopP = 0, pertpopconjsurvived = 0,
                                          timepertpopconjextinct = 0,
                                          timefinal = 1, tmaxshort = 0,
                                          eqreached = 1, infgrowth = 0)
                  }
                } else {
                  # No simulations over time performed, so set values to NA
                  abunfinalconj <- list(R = rep(NA, nspecies), Rtotal = NA,
                                        npopR = NA,
                                        P = rep(NA, nspecies), Ptotal = NA,
                                        npopP = NA, pertpopconjsurvived = NA,
                                        timepertpopconjextinct = NA,
                                        timefinal = NA, tmaxshort = NA,
                                        eqreached = NA, infgrowth = NA)
                }
                
                # Total (cells / mL) and relative (fractions of total) abundances
                abuntotalconj <- abunfinalconj$Rtotal + abunfinalconj$Ptotal
                relabunRconjsp[seq_len(nspecies)] <- abunfinalconj$R / abuntotalconj
                relabunPconjsp[seq_len(nspecies)] <- abunfinalconj$P / abuntotalconj
                
                # Using abunfinalconj$R + abunfinalconj$P > smallstate if the
                # abundances are NA leads to 0 instead of NA, so instead use
                # if-else construct.
                if(abunfinalconj$eqreached == 0 | simulateinvasion == FALSE) {
                  nspeciesconj <- NA 
                  fracPformedbypertpop <- NA
                } else {
                  nspeciesconj <- length(which(
                    abunfinalconj$R + abunfinalconj$P > smallstate
                  ))
                  if(abunfinalconj$Ptotal > 0) {
                    fracPformedbypertpop <- sum(abunfinalconj$P[pertpopconj]) /
                      abunfinalconj$Ptotal
                  } else {
                    fracPformedbypertpop <- NA
                  }
                }
                
                indexdatanew <- indexdata + nspecies
                
                data[indexdata:(indexdatanew - 1), ] <- cbind(
                  nsimul, nspecies, abunmodelcode, intmean, selfintmean,
                  cost, conjratecode, taxmatcode, iter, seq_len(nspecies),
                  c(0, abundance),
                  diag(intmat), c(growthrate), newgrowthratecode, iterintmat,
                  matrix(rep(eqinfo, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(eqinfoconj, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(eqinfoecol, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(eqinfoepi, nspecies), nrow = nspecies, byrow = TRUE),
                  compstability,
                  abunfinal$infgrowth, abunfinalconj$infgrowth,
                  abunfinal$eqreached, abunfinalconj$eqreached,
                  abunfinal$tmaxshort, abunfinalconj$tmaxshort,
                  abunfinal$timefinal, abunfinalconj$timefinal,
                  abunfinal$Rtotal, abunfinalconj$Rtotal, abunfinalconj$Ptotal,
                  abuntotalconj, abunfinalconj$Rtotal/abuntotalconj,
                  abunfinal$npopR, abunfinalconj$npopR, abunfinalconj$npopP,
                  abunfinalconj$npopR + abunfinalconj$npopP,
                  abunfinalconj$Ptotal / abuntotalconj,
                  nspeciesconj, abunfinalconj$npopR / nspeciesconj,
                  abunfinalconj$npopP / nspeciesconj,
                  fracPformedbypertpop, abunfinalconj$pertpopconjsurvived,
                  abunfinalconj$timepertpopconjextinct,
                  matrix(rep(relabunRsp, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(relabunRconjsp, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(relabunPconjsp, nspecies), nrow = nspecies, byrow = TRUE),
                  matrix(rep(relabunRconjsp + relabunPconjsp, nspecies),
                         nrow = nspecies, byrow = TRUE)
                )
                indexdata <- indexdatanew
              }
            }
          }
        }
        datatotal[indexdatatotal:(indexdatatotal + nrowdata - 1), ] <- data
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
                            "compstability",
                            "infgrowth", "infgrowthconj",
                            "eqreached", "eqreachedconj",
                            "tmaxshort", "tmaxshortconj",
                            "timefinal", "timefinalconj",
                            "abuntotal", "abunRtotalconj", "abunPtotalconj",
                            "abuntotalconj", "relabunRconj",
                            "npopR", "npopRconj", "npopPconj", "npopconj",
                            "relabunPconj", "nspeciesconj",
                            "fracspeciesRconj", "fracspeciesPconj",
                            "fracPformedbypertpop", "pertpopconjsurvived",
                            "timepertpopconjextinct",
                            paste0("relabunRsp", seq_len(maxnspecies)),
                            paste0("relabunRconjsp", seq_len(maxnspecies)),
                            paste0("relabunPconjsp", seq_len(maxnspecies)),
                            paste0("relabunconjsp", seq_len(maxnspecies)))
        
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
        
        if(simulateinvasion == TRUE) {
          # If invasion was simulated, get summary of generated data for all
          # combinations of costs, conjugation rates, and taxmatcodes
          summarydatasimulation <- as_tibble(data) %>%
            group_by(cost, conjratecode, taxmatcode) %>%
            summarise(
              across(c(infgrowth, infgrowthconj, eqreached, eqreachedconj,
                       tmaxshort, tmaxshortconj, pertpopconjsurvived),
                     getfracnotzero, .names = "{.col}{.fn}"),
              timefinalmedian = median(timefinal),
              timefinalconjmedian = median(timefinalconj),
              across(c(abuntotal, abunRtotalconj,
                       starts_with("relabunR"), abunPtotalconj,
                       starts_with("relabunP"), abuntotalconj,
                       relabunRconj, starts_with("relabunconjsp"),
                       starts_with("npop"), nspeciesconj,
                       fracspeciesRconj, fracspeciesPconj, fracPformedbypertpop,
                       timepertpopconjextinct),
                     getsummary4, .names = "{.col}{.fn}"),
              .groups = "drop"
            )
          summarydata <- full_join(summarydata, summarydatasimulation,
                                   by = c("cost", "conjratecode", "taxmatcode"))
        }
        
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
colnames(datatotal) <- colnames(data)

summary(warnings())
if(simulateinvasion == TRUE) {
  eqnotreached <- 1 - plotdata[, "eqreachedfrac"]
  eqnotreachedconj <- 1 - plotdata[, "eqreachedconjfrac"]
  if(any(eqnotreached > 0)) {
    warning("Fraction of simulations where equilibrium has not been reached",
            " after pertubation with plasmid-free\nbacteria ranges from ",
            signif(summary(eqnotreached)["Min."], 3), " to ",
            signif(summary(eqnotreached)["Max."], 3)," (mean: ",
            signif(summary(eqnotreached)["Mean"], 3), "; median: ",
            signif(summary(eqnotreached)["Median"], 3), ").",
            " Use silenteqnotreached = FALSE in perturbequilibrium_pinnew()",
            " for more info")
    plot(density(eqnotreached), main = "Density of not reaching equilibria")
    grid()
  }
  if(any(eqnotreachedconj > 0)) {
    warning("Fraction of simulations where equilibrium has not been reached ",
            " after pertubation with plasmid-bearing\nbacteria ranges from ",
            signif(summary(eqnotreachedconj)["Min."], 3), " to ",
            signif(summary(eqnotreachedconj)["Max."], 3), " (mean: ",
            signif(summary(eqnotreachedconj)["Mean"], 3), "; median: ",
            signif(summary(eqnotreachedconj)["Median"], 3), "). ",
            " Use silenteqnotreached = FALSE in perturbequilibrium_pinnew()",
            " for more info")
    plot(density(eqnotreachedconj),
         main = "Density of not reaching\nequilibria with conjugation")
    grid()
  }
}


#### Saving settings and output to CSV and RDS files ####
DateTimeStamp <- paste0(format(Sys.time(), format = "%Y_%m_%d_%H_%M"),
                        "_PInNewSp")
if(PReplMostAbun == FALSE) {
  DateTimeStamp <- paste0(DateTimeStamp, "_PReplLeastAbun")
}

# Saving the data as R-object into an R data file takes much less space than
# saving it as csv. On reading the data back into R, see the section 'Reading
# previously saved data from a CSV file' below.
saveRDS(object = plotdata, file = paste0(DateTimeStamp, "_plotdata.rds"))
saveRDS(object = datatotal, file = paste0(DateTimeStamp, "_datatotal.rds"))

if(nrow(plotdata) > 250000) {
  warning("Not saved 'plotdata' to CSV-file because the number of rows (",
          nrow(plotdata), ") exceeds 250000.")
} else {
  write.csv(plotdata, file = paste0(DateTimeStamp, "_plotdata.csv"),
            quote = FALSE, row.names = FALSE)
}
names(conjrateset) <- paste0("conjrateset", seq_along(conjrateset))
settings <- c(list(nsimul = nsimul, niterintmat = niterintmat,
                   simulateinvasion = simulateinvasion,
                   smallstate = smallstate, smallchange = smallchange,
                   tstep = formals(perturbequilibrium_pinnew)$tstep,
                   saveplots = saveplots, nspeciesset = nspeciesset,
                   abunmodelset = abunmodelset, totalabun = totalabun,
                   intmeanset = intmeanset, selfintmeanset = selfintmeanset,
                   newgrowthratecode = newgrowthratecode,
                   costset = costset, conjrateset, taxmattype = taxmattypeset,
                   PFrom = "PInNewSp", PReplMostAbun = PReplMostAbun,
                   duration = duration))
for(index in seq_along(settings)) {
  # Using write.table instead of write.csv() to be able to use append = TRUE
  write.table(t(as.data.frame(settings[index])), 
              file = paste0(DateTimeStamp, "_settings.csv"), append = TRUE,
              quote = FALSE, sep = ",", col.names = FALSE)
}
capture.output(sessionInfo(),
               file = paste0(DateTimeStamp, "_sessioninfo_base.txt"))
if(requireNamespace("sessioninfo")) {
  capture.output(sessioninfo::session_info(),
                 file = paste0(DateTimeStamp, "_sessioninfo.txt"))
}


#### Reading previously saved data from R data files or from a CSV file ####
# Put the R data files in the working directory:
# getwd()                      # Prints the current working directory
# list.files(pattern = "data") # Prints file names that contain 'data'.

# Then provide the DateTimeStamp corresponding to the desired file (which can be
# copied from the filenames printed above):
# DateTimeStamp <- "YYYY_mm_dd_HH_MM_PInNewSp"

# Then read the R data files into R by choosing the appropriate section below.

##### Reading from R data files #####
# plotdata <- as.data.frame(
#   readRDS(file = file.path(getwd(), paste0(DateTimeStamp, "_plotdata.rds"))))
# datatotal <- as.data.frame(
#   readRDS(file = file.path(getwd(), paste0(DateTimeStamp, "_datatotal.rds"))))
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))

# The error 'cannot open compressed file 'YYYY_mm_dd_HH_MM_plotdata.rds',
# probable reason 'No such file or directory' occurs if the file is not in the
# working directory (which can be checked by using list.files()) or the
# extension '.rds' is not present in the filename of the actual file. In the
# latter case, either remove '.rds' from the filename supplied to readRDS() or
# add '.rds' to the filename of the actual file.

##### Reading from csv-files #####
# plotdata <- read.csv(
#   file = file.path(getwd(), paste0(DateTimeStamp, "_plotdata.csv")),
#   header = TRUE, sep = ",", quote = "\"", dec = ".", stringsAsFactors = FALSE)
# If plotdata has only one column, probably a semicolon instead of a comma was
# used as separator in the CSV file. So read the file again using that separator.
# if(near(ncol(plotdata), 1L)) {
#   plotdata <- read.csv(
#     file = file.path(getwd(), paste0(DateTimeStamp, "_plotdata.csv")),
#     header = TRUE, sep = ";", quote = "\"", dec = ".", stringsAsFactors = FALSE)
# }
# nspeciesset <- sort(unique(plotdata[, "nspecies"]))


#### Labels and limits for plots ####
labspecies <- paste("Sp.", seq_len(maxnspecies))
names(labspecies) <- seq_len(maxnspecies)
labnspecies <- paste(nspeciesset, "species")
names(labnspecies) <- nspeciesset
# abunmodel_lab has been defined and named in the settings before the simulations.
labcost <- paste0("Fitness cost\n", costset, "/h")
names(labcost) <- costset
labconjrate <- paste("Conjset", seq_along(conjrateset))
names(labconjrate) <- seq_along(conjrateset)
labtaxmat <- c("all\nconjugation\nrates equal",
               "lower\nintersp. conj.\nrates initP")
names(labtaxmat) <- seq_along(taxmattypeset)
# '.multi_line = FALSE' to collapse facet labels into a single label
mylabeller <- labeller(species = labspecies, nspecies = labnspecies,
                       abunmodelcode = abunmodel_lab,
                       cost = labcost, conjratecode = labconjrate,
                       taxmatcode = labtaxmat, .multi_line = FALSE,
                       .default = label_value)
memory.limit(size = 10000)
plotdata <- as.data.frame(plotdata)
datatotal <- as.data.frame(datatotal)
limitsfraction <- c(0, 1)
newgrowthratecode <- unique(plotdata$newgrowthratecode)
if(length(newgrowthratecode) > 1L) {
  stop("'newgrowthratecode' should have length 1!")
}
# Round the limits to one decimal place, while ensuring that all the data is
# within the rounded limits.
limitsgrowthrate <- c(floor(min(plotdata[, "growthratemin"])*10)/10,
                      ceiling(max(plotdata[, "growthratemax"])*10)/10)
stat_type <- c("min", "mean", "median", "max")
names(stat_type) <- c("Min.", "Mean", "Median", "Max.")


#### Plot output ####
# If the error '$ operator is invalid for atomic vectors' arises, the matrix
# 'plotdata' has not yet been converted to a dataframe, run the next line to do
# so:
# plotdata <- as.data.frame(plotdata)


## Compare equilibrium characteristics for the models without and with plasmids.
# If invasion has been simulated, data on infinite growth and reaching
# equilibrium after perturbation is also plotted.
CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction,
           facetx = "abunmodelcode", facety = "nspecies",
           filename = "fracstablefewfacets")
CreatePlot(fillvar = "fracstable", filltitle = "Fraction stable",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(fillvar = "fracstableconj",
           filltitle = "Fraction stable\nwith conjugation",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "1 - fracstable", filltitle = "Fraction unstable",
           filltype = "continuous", limits = limitsfraction,
           facetx = "abunmodelcode", facety = "nspecies",
           filename = "fracunstablefewfacets")
CreatePlot(fillvar = "1 - fracstable", filltitle = "Fraction unstable",
           filltype = "continuous", limits = limitsfraction,
           filename = "fracunstablecontinuous")
CreatePlot(fillvar = "1 - fracstableconj",
           filltitle = "Fraction unstable\nwith conjugation",
           filltype = "continuous", limits = limitsfraction,
           filename = "fracunstableconjcontinuous")

CreatePlot(fillvar = "fracstableecol",
           filltitle = "Fraction of simulations that\nis ecologically stable",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "fracstableepi",
           filltitle = "Fraction of simulations that\nis epidemiologically\nstable",
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "1 - fracstableecol",
           filltitle = "Fraction of simulations that\nis ecologically unstable",
           filltype = "continuous", limits = limitsfraction)
CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "1 - fracstableecol",
           filltitle = "Fraction of simulations that\nis ecologically unstable",
           filltype = "continuous", limits = limitsfraction, tag = "B",
           palette = "turbo", filename = "FigS03B", width = 1.15 * 1127,
           height = 1.1 * 2756)
CreatePlot(fillvar = "1 - fracstableepi",
           filltitle = "Fraction of simulations that\nis epidemiologically unstable",
           filltype = "continuous", limits = limitsfraction)

# Show the effect of adding conjugation on stability
CreatePlot(fillvar = "fracunstableunstable + fracneutralneutral + fracstablestable",
           filltitle = "Fraction stability\nunchanged\nby conjugation",
           filltype = "continuous", limits = limitsfraction)

eq_states <- c("unstable", "neutral", "stable")
for(eq_status_without in eq_states) {
  for(eq_status_with in eq_states) {
    print(CreatePlot(fillvar = paste0("frac", eq_status_without, eq_status_with),
                     filltitle = paste0("Fraction ",
                                        eq_status_without, " without and\n",
                                        eq_status_with, " with conjugation"),
                     filltype = "continuous", limits = limitsfraction))
  }
}

if(simulateinvasion == TRUE) {
  CreatePlot(fillvar = "infgrowthfrac",
             filltitle = "Fraction of simulations where\ninfinite growth occurred",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedfrac",
             filltitle = "Fraction of simulations where\nthe equilibrium was reached",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "tmaxshortfrac",
             filltitle = "Fraction of simulations where\ntmax was too short",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "infgrowthconjfrac",
             filltitle = "Fraction of simulations where\ninfinite growth occurred",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "infgrowthconjfrac",
             filltitle = "Fraction of simulations where\ninfinite growth occurred",
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "turbo", filename = "FigS10B")
  
  CreatePlot(fillvar = "eqreachedconjfrac",
             filltitle = "Fraction of simulations where\nthe equilibrium was reached ",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "eqreachedconjfrac",
             filltitle = "Fraction of simulations where\nthe equilibrium was reached ",
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "turbo", filename = "FigS09B")
  CreatePlot(fillvar = "tmaxshortconjfrac",
             filltitle = "Fraction of simulations where\ntmax was too short",
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "tmaxshortconjfrac",
             filltitle = "Fraction of simulations where\ntmax was too short",
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "turbo", filename = "FigS11B")
}

## Growth rates
# Plot summary data for the calculated growth rates 
for(ind_stat_type in seq_along(stat_type)) {
  print(CreatePlot(fillvar = paste0("growthrate", stat_type[ind_stat_type]),
                   filltitle = paste(names(stat_type[ind_stat_type]), "growth rate"),
                   filltype = "continuous", limits = limitsgrowthrate))
}

pS02 <- CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                             near(conjratecode, 1), near(taxmatcode, 1)),
           fillvar = "growthratemean", filltitle = "Mean growth rate",
           filltype = "continuous", limits = NULL, tag = "B", palette = "mako",
           save = FALSE)
pS02 +
  labs(y = "Mean intraspecies\ninteraction coefficient", tag = NULL) +
  facet_wrap(facets = "nspecies", nrow = 1, labeller = label_both)
ggsave(paste0(DateTimeStamp, "FigS02_growthratecode", newgrowthratecode, ".png"),
       width = 2028, height = 1.1 * 704, dpi = 300, units = "px")

# Show the relation of interactions and species-specific growth rate required to
# obtain an equilibrium. Costs and conjugation rate do not affect growth rate,
# so data has been filtered to have only one value for them.
datatotalfiltercostconj <- filter(datatotal, near(cost, costset[1]),
                                  near(conjratecode, 1), near(taxmatcode, 1))
datatotalfiltercostconjsp1 <- filter(datatotalfiltercostconj,
                                     near(species, 1))
ggplot(data = datatotalfiltercostconjsp1,
       aes(x = intmean, y = growthrate)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom",
        panel.border = element_blank(),
        panel.spacing = unit(3, "pt"),
        strip.background = element_rect(color = NA)) +
  geom_point(size = 0.1) +
  facet_grid(rows = vars(selfintmean), cols = vars(nspecies),
             labeller = mylabeller, drop = TRUE, scales = "fixed") +
  scale_color_viridis_c() +
  labs(title = "Growth rate of species 1",
       caption = paste(nsimul, "simulations")) +
  guides(color = guide_colourbar())
if(saveplots == TRUE) {
  ggsave(paste0(DateTimeStamp, "growthratevsintmean_sp1.png"),
         width = 3300, height = 2675, units = "px", dpi = 300)
}

datatotalfiltercostconj_iter1 <- datatotalfiltercostconj[
  which(near(1, datatotalfiltercostconj[, "iter"])), ]

## Plot summary data on the number of simulations in creating intmat needed to
# find a stable equilibrium with the model without plasmids
for(ind_stat_type in seq_along(stat_type)) {
  print(CreatePlot(fillvar = paste0("iterintmat", stat_type[ind_stat_type]),
                   filltitle = paste(names(stat_type[ind_stat_type]),
                                     "number of\nsimulations to reach\nstable",
                                     "equilibrium"),
                   filltype = "continuous", limits = c(1, niterintmat)))
}

if(simulateinvasion == TRUE) {
  subplasmidfree <- "Perturbation with plasmid-free bacteria"
  subplasmidbearing <- "Perturbation with plasmid-bearing bacteria"
  if(PReplMostAbun == FALSE) {
    title_add <- " (PReplLeastAbun)"
  } else {
    title_add <- NULL
  }
  
  limitstime <- range(plotdata[, "timefinalmedian"],
                      plotdata[, "timefinalconjmedian"], na.rm = TRUE)
  title <- paste0("Time to reach equilibrium after perturbation", title_add)
  CreatePlot(fillvar = "timefinalmedian", filltitle = "Median time",
             filltype = "continuous", limits = limitstime,
             title = title, subtitle = subplasmidfree, rotate_legend = TRUE)
  CreatePlot(fillvar = "timefinalconjmedian", filltitle = "Median time",
             filltype = "continuous", limits = limitstime,
             title = title, subtitle = subplasmidbearing, rotate_legend = TRUE)
  
  ## Plots on survival and extinction after perturbation
  limitsmeannspecies <- range(plotdata[, "npopRmean"],
                              plotdata[, "nspeciesconjmean"], finite = TRUE)
  title <- paste0("Number of species surviving after perturbation", title_add)
  # Note: In the model without plasmids, the number of populations is equal to
  # the number of species.
  CreatePlot(fillvar = "npopRmean", filltitle = "Mean total number\nof species",
             filltype = "continuous", limits = limitsmeannspecies,
             title = title, subtitle = subplasmidfree,
             filename = "nspeciesRmean")
  CreatePlot(fillvar = "nspeciesconjmean", filltitle = "Mean total number\nof species",
             filltype = "continuous", limits = limitsmeannspecies,
             title = title, subtitle = subplasmidbearing)
  
  title <- paste0("Number of species extinct after perturbation", title_add)
  CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                               near(conjratecode, 1), near(taxmatcode, 1)),
             fillvar = "nspecies - npopRmean",
             filltitle = "Mean number of species\nthat went extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidfree,
             facetx = "abunmodelcode", facety = "nspecies",
             filename = "nspeciesRmeanextinctfewfacets")
  CreatePlot(fillvar = "nspecies - npopRmean",
             filltitle = "Mean number of species\nthat went extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidfree,
             filename = "nspeciesRmeanextinct")
  CreatePlot(fillvar = "nspecies - nspeciesconjmean",
             filltitle = "Mean number of species\nthat went extinct",
             filltype = "continuous", limits = c(0, maxnspecies),
             title = title, subtitle = subplasmidbearing,
             filename = "nspeciesconjmeanextinct")
  
  title <- paste0("Fraction of species extinct after perturbation", title_add)
  CreatePlot(dataplot = filter(plotdata, near(cost, costset[1]),
                               near(conjratecode, 1), near(taxmatcode, 1)),
             fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of species\nthat went extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree,
             facetx = "abunmodelcode", facety = "nspecies",
             filename = "fracspeciesextinctmeanfewfacets")
  CreatePlot(fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of species\nthat went extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree,
             filename = "fracspeciesextinctmean")
  CreatePlot(fillvar = "1 - (npopRmean / nspecies)",
             filltitle = "Mean fraction of species\nthat went extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidfree)
  CreatePlot(fillvar = "1 - (nspeciesconjmean / nspecies)",
             filltitle = "Mean fraction of species\nthat went extinct",
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing,
             filename = "fracspeciesextinctmeanconj")
  CreatePlot(fillvar = "1 - (nspeciesconjmean / nspecies)",
             filltitle = "Mean fraction of species\nthat went extinct",
             filltype = "continuous", limits = limitsfraction, tag = "A",
             palette = "plasma", filename = "FigS05A")
  
  title <- paste0("Fraction of species after perturbation", title_add)
  CreatePlot(fillvar = "fracspeciesRconjmean",
             filltitle = paste("Mean fraction of surviving species\nwith a",
                               "plasmid-free population "),
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "npopRconjmean / nspecies",
             filltitle = paste("Mean fraction of initial species\nwith a",
                               "plasmid-free population "),
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "fracspeciesPconjmean",
             filltitle = paste("Mean fraction of surviving species\nwith a",
                               "plasmid-bearing population "),
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "fracspeciesPconjmean",
             filltitle = paste("Mean fraction of surviving species\nwith a",
                               "plasmid-bearing population "),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "plasma", filename = switch(newgrowthratecode,
                                                   "Fig05B", "Fig04B", "Fig06B"))
  
  CreatePlot(fillvar = "npopPconjmean / nspecies",
             filltitle = paste("Mean fraction of initial species\nwith a",
                               "plasmid-bearing population "),
             filltype = "continuous", limits = limitsfraction,
             title = title, subtitle = subplasmidbearing)
  CreatePlot(fillvar = "npopPconjmean / nspecies",
             filltitle = paste("Mean fraction of initial species\nwith a",
                               "plasmid-bearing population "),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "plasma", filename = paste0("FigS05B_growthratecode",
                                                   newgrowthratecode))
  
  ## Plots of fractions of bacteria that are plasmid-free or plasmid-bearing
  # after perturbations. Only abundances where equilibrium was reached are
  # considered.
  title <- paste0("Fraction bacteria after perturbation", title_add)
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("relabunRconj", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of bacteria\nthat is plasmid-free"),
                     filltype = "continuous", limits = limitsfraction,
                     title = title, subtitle = subplasmidbearing))
    
    print(CreatePlot(fillvar =  paste0("relabunPconj", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of bacteria\nthat is plasmid-bearing"),
                     filltype = "continuous", limits = limitsfraction,
                     title = title, subtitle = subplasmidbearing))
    
    print(CreatePlot(fillvar =  paste0("1 - fracPformedbypertpop", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of plasmid-bearing\nbacteria",
                                       "belonging to the\ninitially plasmid-free species"),
                     filltype = "continuous", limits = limitsfraction))
    
    print(CreatePlot(fillvar =  paste0("fracPformedbypertpop", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "fraction of plasmid-bearing\nbacteria",
                                       "belonging to the\ninitially plasmid-bearing species"),
                     filltype = "continuous", limits = limitsfraction))
    
    print(CreatePlot(fillvar =  paste0("timepertpopconjextinct", stat_type[ind_stat_type]),
                     filltitle = paste(names(stat_type[ind_stat_type]),
                                       "time initially plasmid-\nbearing",
                                       "population went\nextinct"),
                     filltype = "continuous", rotate_legend = TRUE))
  }
  
  if(newgrowthratecode == 2) {
    CreatePlot(fillvar =  "relabunPconjmean",
               filltitle = "Mean fraction of bacteria\nthat is plasmid-bearing",
               filltype = "continuous", limits = limitsfraction, tag = "B",
               filename = "Fig03B")
  } else {
    CreatePlot(fillvar =  "relabunPconjmean",
               filltitle = "Mean fraction of bacteria\nthat is plasmid-bearing",
               filltype = "continuous", limits = limitsfraction, tag = "A",
               filename = if(newgrowthratecode == 1) {"Fig05A"} else {"Fig06A"})
  }
  
  CreatePlot(fillvar =  "fracPformedbypertpopmean",
             filltitle = paste("Mean fraction of plasmid-bearing\nbacteria",
                               "belonging to the initially\nplasmid-bearing species"),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             height = 1.05 * 1680, filename = "FigS13B")
  
  CreatePlot(fillvar = "pertpopconjsurvivedfrac",
             filltitle = paste("Fraction of simulations where\nthe initially",
                               "plasmid-bearing\npopulation survived"),
             filltype = "continuous", limits = limitsfraction)
  CreatePlot(fillvar = "pertpopconjsurvivedfrac",
             filltitle = paste("Fraction of simulations where\nthe initially",
                               "plasmid-bearing\npopulation survived"),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "turbo", height = 1.05 * 1680, filename = "FigS12B")
  
  ## Plot total abundances after perturbation with plasmid-free bacteria in
  # models without plasmids. Only abundances where equilibrium was reached are
  # considered. Although costs and conjugation rates do not influence the
  # outcome, they are included as facets to facilitate comparison with plots of
  # abundances after perturbation with plasmid-bearing bacteria.
  
  # log10(0) leads to -Inf and is displayed as a gray square in the plots,
  # leading to different numbers of gray squares between plots. It is also
  # confusing because gray squares can also indicate no data is available
  # because all simulations were discarded because equilibrium was not reached.
  # Therefore I plot log10(1 + x) instead of log10(x). I do not use the built-in
  # function log1p(x) because that uses natural logarithms. For
  # relative abundances I used log(1e-6 + x) because adding 1 is too much.
  title <- paste0("Total abundance after perturbation", title_add)
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abuntotal", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "abundance of\nplasmid-free bacteria)"),
                     filltype = "continuous", title = title,
                     rotate_legend = TRUE, subtitle = subplasmidfree))
  }
  
  ## Plot of total abundances after perturbation with plasmid-bearing bacteria in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abuntotalconj", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "\ntotal abundance)"),
                     filltype = "continuous", title = title, subtitle = subplasmidbearing))
  }
  
  ## Plot total abundances of plasmid-free populations after perturbations in
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abunRtotalconj", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "abundance of\nplasmid-free bacteria)"),
                     filltype = "continuous", title = title, subtitle = subplasmidbearing))
  }
  
  ## Plot total abundances of plasmid-bearing populations after perturbations for
  # models with plasmids. Only abundances where equilibrium was reached are
  # considered.
  for(ind_stat_type in seq_along(stat_type)) {
    print(CreatePlot(fillvar = paste0("log10(1 + abunPtotalconj", stat_type[ind_stat_type], ")"),
                     filltitle = paste("Log10(1 +", names(stat_type[ind_stat_type]),
                                       "abundance of\nplasmid-bearing bacteria)"),
                     filltype = "continuous", title = title,
                     subtitle = subplasmidbearing, rotate_legend = TRUE))
  }
  
  ## Plots comparing species abundances after perturbation with plasmid-bearing
  # bacteria
  limits_mean <- range(c(plotdata[, "relabunRsp1mean"],
                         plotdata[, "relabunconjsp1mean"]), na.rm = TRUE)
  limits_median <- range(c(plotdata[, "relabunRsp1median"],
                           plotdata[, "relabunconjsp1median"]), na.rm = TRUE)
  filltitle_mean <- paste("Mean rel. abundance\nof sp1 after perturbation",
                          "with R\nof newly introduced species 1")
  filltitle_conjmean <- paste("Mean rel. abundance\nof sp1 after perturbation",
                              "with P\nof newly\nintroduced species 1")
  filltitle_P1conjmean <- paste("Mean rel. abundance\nof P1 after perturbation",
                                "with P\nof newly introduced species 1")
  filltitle_median <- paste("Median rel. abundance\nof sp1 after perturbation",
                            "with R\nof newly introduced species 1")
  filltitle_conjmedian <- paste("Median rel. abundance\nof sp1 after perturbation",
                                "with P\nof newly introduced species 1")
  
  # Note: relabunRspXmean is the mean of the total population because no
  # plasmids are present (and thus no plasmid-bearing population).
  CreatePlot(fillvar = "relabunRsp1mean", filltitle = filltitle_mean,
             filltype = "continuous", limits = limitsfraction, rotate_legend = TRUE,
             filename = "relabunRsp1meancontinuous_alt")
  CreatePlot(fillvar = "relabunRsp1mean", filltitle = filltitle_mean,
             filltype = "continuous", rotate_legend = TRUE,
             filename = "relabunRsp1meancontinuousnolim")
  CreatePlot(fillvar = "relabunRsp1mean", filltitle = filltitle_mean,
             filltype = "continuous", limits = limits_mean, rotate_legend = TRUE,
             filename = "relabunRsp1meancontinuouschangedlim")
  CreatePlot(fillvar = "log10(1e-6 + relabunRsp1mean)",
             filltitle = paste0("Log10(1e-6 + ", filltitle_mean, ")"),
             filltype = "continuous",
             filename = "relabunRsp1meancontinuousloglim")
  
  CreatePlot(fillvar = "relabunconjsp1mean", filltitle = filltitle_conjmean,
             filltype = "continuous", limits = limitsfraction, rotate_legend = TRUE,
             filename = "relabunconjsp1meancontinuous_alt")
  CreatePlot(fillvar = "relabunconjsp1mean", filltitle = filltitle_conjmean,
             filltype = "continuous", rotate_legend = TRUE,
             filename = "relabunconjsp1meancontinuousnolim")
  CreatePlot(fillvar = "relabunconjsp1mean", filltitle = filltitle_conjmean,
             filltype = "continuous", limits = limits_mean, rotate_legend = TRUE,
             filename = "relabunconjsp1meancontinuouschangedlim")
  CreatePlot(fillvar = "log10(1e-6 + relabunconjsp1mean)",
             filltitle = paste0("Log10(1e-6 + ", filltitle_conjmean, ")"),
             filltype = "continuous",
             filename = "relabunconjsp1meancontinuousloglim")
  
  CreatePlot(fillvar = "relabunPconjsp1mean", filltitle = filltitle_P1conjmean,
             filltype = "continuous", limits = limitsfraction, rotate_legend = TRUE,
             filename = "relabunPconjsp1meancontinuous")
  CreatePlot(fillvar = "relabunPconjsp1mean", filltitle = filltitle_P1conjmean,
             filltype = "continuous", rotate_legend = TRUE,
             filename = "relabunPconjsp1meancontinuousnolim")
  CreatePlot(fillvar = "relabunPconjsp1mean", filltitle = filltitle_P1conjmean,
             filltype = "continuous", limits = limits_mean, rotate_legend = TRUE,
             filename = "relabunPconjsp1meancontinuouschangedlim")
  CreatePlot(fillvar = "log10(1e-6 + relabunPconjsp1mean)",
             filltitle = paste0("Log10(1e-6 + ", filltitle_P1conjmean, ")"),
             filltype = "continuous",
             filename = "relabunPconjsp1meancontinuousloglim")
  
  CreatePlot(fillvar = "relabunRsp1median", filltitle = filltitle_median,
             filltype = "continuous", limits = limitsfraction,
             filename = "relabunRsp1mediancontinuouschangedlim")
  CreatePlot(fillvar = "relabunconjsp1median", filltitle = filltitle_conjmedian,
             filltype = "continuous", limits = limitsfraction,
             filename = "relabunconjsp1mediancontinuouschangedlim")
  
  if(PReplMostAbun == TRUE) {
    add_filltitle <- paste0("after\nreplacing some R of the most-abundant",
                            " sp.\nwith R of the newly introduced sp.")
    add_filltitleconj <- paste0("after\nreplacing some R of the most-abundant",
                                " sp.\nwith P of the newly introduced sp.")
  } else {
    add_filltitle <- paste0("after\nreplacing some R of the least-abundant",
                            " sp.\nwith R of the newly introduced sp.")
    add_filltitleconj <- paste0("after\nreplacing some R of the least-abundant",
                                " sp.\nwith P of the newly introduced sp.")
  }
  for(species_i in unique(c(1, nspeciesset)))  {
    for(ind_stat_type in seq_along(stat_type)) {
      print(CreatePlot(fillvar = paste0("relabunRsp", species_i,
                                        stat_type[ind_stat_type]),
                       filltitle = paste0(names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitle),
                       filltype = "continuous", limits = limitsfraction))
      
      print(CreatePlot(fillvar = paste0("log10(1e-6 + relabunRsp", species_i,
                                        stat_type[ind_stat_type], ")"),
                       filltitle = paste0("Log10(1e-6 + ",
                                          names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitle, ")"),
                       filltype = "continuous", rotate_legend = TRUE))
      
      print(CreatePlot(fillvar = paste0("relabunconjsp", species_i,
                                        stat_type[ind_stat_type]),
                       filltitle = paste0(names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitleconj),
                       filltype = "continuous", limits = limitsfraction))
      
      print(CreatePlot(fillvar = paste0("log10(1e-6 + relabunconjsp", species_i,
                                        stat_type[ind_stat_type], ")"),
                       filltitle = paste0("Log10(1e-6 + ",
                                          names(stat_type[ind_stat_type]),
                                          " rel. abundance of sp", species_i, " ",
                                          add_filltitleconj, ")"),
                       filltype = "continuous", rotate_legend = TRUE))
    }
  }
}

if(PReplMostAbun && simulateinvasion) {
if(newgrowthratecode == 1) {
  CreatePlot(fillvar = "relabunconjsp1mean",
             filltitle = paste0("Mean relative abundance of the\ninitially",
                                " plasmid-bearing species "),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "rocket", filename = "FigS14B_growthrate1")
}

if(newgrowthratecode == 2) {
  CreatePlot(fillvar = "relabunconjsp1mean",
             filltitle = paste0("Mean relative abundance of the\ninitially",
                                " plasmid-bearing species "),
             filltype = "continuous", limits = limitsfraction, tag = "A",
             palette = "rocket", height = 1.05 * 1680, filename = "FigS13A")
}

if(newgrowthratecode == 3) {
  CreatePlot(fillvar = "relabunconjsp1mean",
             filltitle = paste0("Mean relative abundance of the\ninitially",
                                " plasmid-bearing species "),
             filltype = "continuous", limits = limitsfraction, tag = "B",
             palette = "rocket", filename = "FigS14B_growthrate3")
}

CreatePlot(fillvar = "1 - relabunconjsp1mean",
           filltitle = paste0("Mean relative abundance of the\ninitially",
                              " plasmid-bearing species ", add_filltitleconj),
           filltype = "continuous", limits = limitsfraction)

CreatePlot(fillvar = "log10(1 - relabunconjsp1mean + 1e-6)",
           filltitle = paste0("log10(1e-6 + Mean relative abundance of\nthe",
                              " initially plasmid-free species) ",
                              add_filltitleconj),
           filltype = "continuous", limits = NULL)
}
