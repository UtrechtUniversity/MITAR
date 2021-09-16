######################################################
## Microbiome analysis in the MITAR/BEWARE project  ##
######################################################

#### Introduction ####
# Following SVEPM workshop by Kers and Fischer, available at
# https://an-jg.github.io/SVEPM2021_Workshop/


#### Setting file path ####
pathseqdata <- "D:/Userdata/Jesse Nieuw/Documents/SurfDriveF/TransmExp/Proef/SequenceData/"

#### Loading libraries ####
library(phyloseq)
# library(microbiome) # read_phyloseq()
library(data.table) # datatable
library(DT) # datatable()

#### Reading data ####
# The sequencing data was converted to a phyloseq object and then saved as .rds-file, using
# R version 4.0.2 (2020-06-22) on a 64-bit x86_64-pc-linux-gnu running under Ubuntu 18.04.5 LTS.
# Now this rds-file is read again to get the phyloseq object

# QUESTION: taxonomic data was also added (from SILVA?), but is NOT YET A TREE ?
# In the workshop a tree-file is added

ps <- readRDS(file = paste0(pathseqdata, "20210203_Egil180x_SILVA_v3v4_merged.rds"))
ps

head(datatable(tax_table(ps))) # I get the raw sequences NOT the taxa
## ERROR: no tree is in the phyloseq object

ps1 <- subset_taxa(ps, Domain!="NA")
summarize_phyloseq(ps)

