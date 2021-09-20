######################################################
## Microbiome analysis in the MITAR/BEWARE project  ##
######################################################


#### Introduction ####
# Analysis of microbiome data from transmission experiments in the MITAR/BEWARE
# project.


#### References ####
# I used code from a SVEPM2021 workshop by Kers and Fischer
# (https://an-jg.github.io/SVEPM2021_Workshop/), from the phyloseq tutorials by
# Paul J. McMurdie (https://joey711.github.io/phyloseq/), and from the DADA2
# tutorial by Benjamin Callahan (http://benjjneb.github.io/dada2/tutorial.html)


#### To do ####
# Check if days indicate the age of chicken, OR number of days after inoculation?

# Based on id_group, add info on treatment to sample_data.


#### Setting file path ####
# Path should point to the folder where the phyloseq oject is saved, not to the file itself
pathseqdata <- "D:/Userdata/Jesse Nieuw/Documents/SurfDriveF/TransmExp/Proef/SequenceData/"
# pathseqdata <- "C:/Users/3501477/surfdrive/TransmExp/Proef/SequenceData/"


#### Loading libraries ####
library(phyloseq)   # read and analyse phyloseq objects
library(microbiome) # summarize_phyloseq(), transform()
# library(data.table) # datatable()
# library(DT)         # datatable()


#### Reading and inspecting phylosec object ####
# The sequencing data was converted to a phyloseq object and then saved as .rds-file, using
# R version 4.0.2 (2020-06-22) on a 64-bit x86_64-pc-linux-gnu running under Ubuntu 18.04.5 LTS.
# Now this rds-file is read again to get the phyloseq object in R.

ps <- readRDS(file = paste0(pathseqdata, "20210203_Egil180x_SILVA_v3v4_merged.rds"))
ps 
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8055 taxa and 177 samples ]
# sample_data() Sample Data:       [ 177 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 8055 taxa by 6 taxonomic ranks ]

# So ps is a phyloseq-object containing an otu_table, a sample_data, and a
# tax_table. There is no phy_tree() object (i.e., the phy_tree slot is empty),
# because the taxonomic data is not yet assembled in a taxonomic tree.


#### Inspecting otu_table ####
# The otu_table consists of 177 rows which names are the 177 sample identifiers,
# and 8055 columns which names are the 8055 identified OTUs. The numbers in the
# table give the number of sequence reads identified as those OTUs in each of
# the samples. Therefor sum(otu_table(ps)) equals the 10262667 reads.
# Below the structure of the otu_table is illustrated using a subset containing
# only the first 20 samples with the first 10 taxa, without the column names,
# because the column names are the actual sequences which are hundreds of
# characters long.

ps_otu_table <- otu_table(ps)
dim(ps_otu_table)
# [1]  177 8055

ps_otu_table_subset <- ps_otu_table[1:20, 1:10]
attributes(ps_otu_table_subset) # showing subset
# $dim
# [1] 20 10
# 
# $taxa_are_rows
# [1] FALSE
# 
# $class
# [1] "otu_table"
# attr(,"package")
# [1] "phyloseq"
# 
# $dimnames
# $dimnames[[1]]
# [1] "A_1-10_G_1_D_14_S103"   "A_1-11_G_1_D_5_S127"    "A_1-12_G_1_D_5_S146"    "A_1-13_G_1_D_14_S134"   "A_1-14_G_1_D_14_S89"    "A_1-15_G_1_D_5_S106"   
# [7] "A_1-2_G_1_D_14_S62"     "A_1-3_G_1_D_14_S60"     "A_1-4_G_1_D_14_S94"     "A_1-5_G_1_D_5_S58"      "A_1-6_G_1_D_5_S135"     "A_1-8_G_1_D_14_S32"    
# [13] "A_10-1_G_10_D_5_S85"    "A_10-10_G_10_D_14_S27"  "A_10-11_G_10_D_5_S68"   "A_10-12_G_10_D_5_S79"   "A_10-13_G_10_D_14_S116" "A_10-14_G_10_D_5_S20"  
# [19] "A_10-15_G_10_D_14_S163" "A_10-2_G_10_D_14_S13"  
# 
# $dimnames[[2]]
# [1] "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGT[...]CGTGGG" 
# [2] "GGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAAC[...]CGTGGG"
# [3] "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGT[...]CGTGGG" 
# [4] "GGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGATGACGGCCTTCGGGTTGTAAAGC[...]CGTGGG"                          
# [5] "GGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAAC[...]CGTGGG" 
# [6] "GGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACT[...]CGTGGG"                       
# [7] "GGAATATTGCACAATGGGGGAAACCCTGATGCAGCGATGCCGCGTGAAGGAAGAAGTATCTCGGTATGTAAACT[...]CGTGGG"                          
# [8] "GGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGCGAAGAAGTATTTCGGTATGTAAAGC[...]CGTGGG"                          
# [9] "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGT[...]CGTGGG" 
# [10]"GGAATATTGCACAATGGGCGAAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGT[...]CGTGGG"   

# The column names are the actual sequences and thus hundreds of characters long
colnames(ps_otu_table)[1]
plot(table(nchar(colnames(ps_otu_table))))

# To get clearer view on the structure of the otu table, remove the column names
# and show the first 10 taxa for the first 20 samples
ps_otu_table_subset_nocolnames <- ps_otu_table_subset
colnames(ps_otu_table_subset_nocolnames) <- NULL
ps_otu_table_subset_nocolnames

# OTU Table:          [10 taxa and 20 samples]
# taxa are columns
#                         [,1]  [,2] [,3]  [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# A_1-10_G_1_D_14_S103    5266  1916 4632     5  174 1717    3    8 1561  1371
# A_1-11_G_1_D_5_S127     6310     7 5525  4776  529    6    6    2 1732  1658
# A_1-12_G_1_D_5_S146     3673     3 3126  3480  375    7    0    6 1082  1000
# A_1-13_G_1_D_14_S134    5380   237 4719     6   70 2909    0    0 1487  1398
# A_1-14_G_1_D_14_S89     3848  1204 3127    12  139 1845    6    5 1009   970
# A_1-15_G_1_D_5_S106     9737     9 8796  4102  949   12    4    5 2722  2485
# A_1-2_G_1_D_14_S62     11796 26353 4615   186  269 2251    2    6 1032   946
# A_1-3_G_1_D_14_S60      5355   717 4518     4  122 3544    1    0 1617  1421
# A_1-4_G_1_D_14_S94      3722  1496 3086    12  270 1711    8    6 1081   986
# A_1-5_G_1_D_5_S58       7895    18 6179 10148 2197   16    3    8 2728  2875
# A_1-6_G_1_D_5_S135      6640     6 6005  3052  208    1    8    8 1784  1688
# A_1-8_G_1_D_14_S32      3918   845 3180     7  315 3880    0    0 1174  1113
# A_10-1_G_10_D_5_S85     2144  6405 1862  2013    2    0    0    5  626   576
# A_10-10_G_10_D_14_S27   4514  7751 3693     8   19  588 3892 2530 1382  1266
# A_10-11_G_10_D_5_S68    3619  6388 3040  2538   57    2    7    2 1061   955
# A_10-12_G_10_D_5_S79    2711 10188 2270  2531   16    4   10   17  842   741
# A_10-13_G_10_D_14_S116  7725  3960 6679    33   17  824 1415  539 2228  2042
# A_10-14_G_10_D_5_S20    4655  9372 4073  3212   39    5    6    7 1293  1296
# A_10-15_G_10_D_14_S163  3793  6711 3104     8   27  424 5318 1519 1041  1031
# A_10-2_G_10_D_14_S13    6002  2554 5125    35    5  472 1948  768 1806  1589


#### Inspecting sample_data ####
# The sample_data is a phyloseq dataframe. It consists of 177 rows which names
# are the 177 sample identifiers, and a single column named "SampleNames" which
# gives the 177 sample identifiers. These identifiers should be split over
# multiple columns to extract the metadata contained in the sample identifier.

ps_sample_data <- sample_data(ps)
dim(ps_sample_data)
# [1] 177   1

ps_sample_data_subset <- ps_sample_data[1:20, ]
attributes(ps_sample_data_subset) # showing subset
# $names
# [1] "SampleNames"
# 
# $row.names
# [1] "A_1-10_G_1_D_14_S103"   "A_1-11_G_1_D_5_S127"    "A_1-12_G_1_D_5_S146"    "A_1-13_G_1_D_14_S134"   "A_1-14_G_1_D_14_S89"    "A_1-15_G_1_D_5_S106"   
# [7] "A_1-2_G_1_D_14_S62"     "A_1-3_G_1_D_14_S60"     "A_1-4_G_1_D_14_S94"     "A_1-5_G_1_D_5_S58"      "A_1-6_G_1_D_5_S135"     "A_1-8_G_1_D_14_S32"    
# [13] "A_10-1_G_10_D_5_S85"    "A_10-10_G_10_D_14_S27"  "A_10-11_G_10_D_5_S68"   "A_10-12_G_10_D_5_S79"   "A_10-13_G_10_D_14_S116" "A_10-14_G_10_D_5_S20"  
# [19] "A_10-15_G_10_D_14_S163" "A_10-2_G_10_D_14_S13"  
# 
# $.S3Class
# [1] "data.frame"
# 
# $class
# [1] "sample_data"
# attr(,"package")
# [1] "phyloseq"

# Concise enough to show complete table
ps_sample_data

## Obtaining data from sample identifiers
# NOTE: HARDCODED column indices. Attempt to use regular expression within
# strsplit() to split when a digit is followed by an underscore did NOT work,
# because the last digits were discarded during the split: strsplit(id_sample[1],
# split = "[0123456789]_", fixed = FALSE) drops digits from id_sample[1].
# The identifiers of the non-spike and spike samples are coded differently, so
# they are processed separately.

# In the sample identifiers of the non-spike samples, animal ID (which also
# contains the group ID) is preceded by A_, the group ID is preceded by G_, the
# day of sample collection is preceded by D_, and the sample number is preceded
# by _S. 

total_id <- unlist(unname(ps_sample_data))
index_spike <- grep("PBS_spike", total_id, value = FALSE, fixed = TRUE)

id_sample <- total_id[-index_spike] # selecting only non-spike samples
id_sample_split <- unname(t(as.data.frame(strsplit(id_sample, split = "_", fixed = TRUE))))
sample_df <- id_sample_split[, c(2, 4, 6, 7)]
colnames(sample_df) <- c("id_animal", "id_group", "day", "sample_nr")

id_spike <- total_id[index_spike] # selecting only spike samples
id_spike_split <- unname(t(as.data.frame(strsplit(id_spike, split = "_", fixed = TRUE))))
id_animal_spike <- paste(id_spike_split[, 2], id_spike_split[, 3], sep = "_")
# Get id_group from id_animal because groups were not coded separately.
# NOTE: assuming single-digit group number.
id_group_spike <- substr(id_animal_spike, start = 1, stop = 7)
day_spike <- rep("spike", length(id_animal_spike))
sample_nr_spike <- id_spike_split[, 4]
spike_df <- cbind(id_animal = id_animal_spike, id_group = id_group_spike,
                  day = day_spike, sample_nr = sample_nr_spike)

# Merge data (same order as original data) and add it to the phyloseq object
total_sample_data_df <- as.data.frame(rbind(sample_df, spike_df),
                                      row.names = rownames(ps_sample_data))
head(total_sample_data_df)
tail(total_sample_data_df)
# Convert dataframe to phyloseq object to merge correctly
ps <- merge_phyloseq(ps, sample_data(total_sample_data_df))
head(sample_data(ps))
tail(sample_data(ps))
ps


#### Inspecting tax_table ####
# The tax_table consists of 8055 rows which names are the 8055 identified OTUs,
# and 6 columns which names are the taxonomic levels from Kingdom to Genus. 
# The table gives each of the taxonomic levels for all sequences. If a sequence
# could not be assigned to a taxonomic level, it is listed as NA. In a separate
# output ALL missing (NA) taxa at level X where filled with one level higher
# taxonomic data (X-1).

# Below the structure of the tax_table is illustrated using a subset containing
# only the first 20 samples with the first 10 taxa, without the column names,
# because the column names are the actual sequences which are hundreds of
# characters long.

ps_tax_table <- tax_table(ps)
dim(ps_tax_table)
# [1] 8055    6

ps_tax_table[1:3, ] # Showing subset
attributes(ps_tax_table[1:3, ]) # Showing subset
# $dim
# [1] 3 6
# 
# $class
# [1] "taxonomyTable"
# attr(,"package")
# [1] "phyloseq"
# 
# $dimnames
# $dimnames[[1]]
# [1] "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTA[...]GCGTGGG" 
# [2] "GGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAACT[...]CGTGGG"
# [3] "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTA[...]CGTGGG" 
# 
# $dimnames[[2]]
# [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"  

ps_tax_table_norownames <- ps_tax_table
rownames(ps_tax_table_norownames) <- NULL
ps_tax_table_norownames[1:3, ] # Showing subset without rownames
# Taxonomy Table:     [3 taxa by 6 taxonomic ranks]:
#     Kingdom    Phylum           Class                 Order              Family               Genus                 
# sp1 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Escherichia/Shigella"
# sp2 "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"  "Enterococcaceae"    "Enterococcus"        
# sp3 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales" "Enterobacteriaceae" "Escherichia/Shigella"


#### Checks on data ####
summarize_phyloseq(ps)

# Checking for NAs in taxonomic data
taxNA <- as.data.frame(matrix(data = NA, nrow = 6, ncol = 3, byrow = FALSE,
                              dimnames = list(NULL, c("TaxRank", "N_NA", "percent_NA"))))
for(index in 1:dim(ps_tax_table)[2]) {
  taxNA[index, "TaxRank"] <- colnames(ps_tax_table[, index])
  taxNA[index, "N_NA"] <- as.numeric(length(which(is.na(ps_tax_table[, index]))))
}
taxNA[, "percent_NA"] <- round(100*taxNA[, "N_NA"] / dim(ps_tax_table)[1], 2)
taxNA  # Gives same numbers as the table in the read-me file
#   TaxRank N_NA percent_NA
# 1 Kingdom    8       0.10
# 2  Phylum  156       1.94
# 3   Class  255       3.17
# 4   Order  711       8.83
# 5  Family 1509      18.73
# 6   Genus 3821      47.44

# Checking for non-prokaryotic sequences in taxonomic data: One sequence
# assigned to an Archaeon (genus Halococcus), and 94 sequences to Eukaryota (all
# taxa apart from Kingdom listed as "NA")
ps_tax_table_norownames[which(ps_tax_table[, "Kingdom"] != "Bacteria"), ]
table(ps_tax_table[, "Kingdom"])
# Archaea  Bacteria Eukaryota 
# 1      7952        94 


#### Plotting data ####

# Bar plots of composition
ps_compo <- transform(ps, "compositional")
ps_compo_filt = filter_taxa(ps_compo, function(x) sum(x) > .01, prune = TRUE)  # filtering
ps_phylum = tax_glom(ps_compo_filt, taxrank = "Phylum", NArm = FALSE)  # Taxonomic level Phylum - Family
plot_bar(ps_phylum, fill = "Phylum") + facet_grid("day" ~ "id_group", scales = "free_x") 

# Facet_grid currently does NOT result in multiple facets (because columns are
# character instead of factor type?

# Richness
plot_richness(ps, "id_group", "day", measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) 
