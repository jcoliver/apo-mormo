# Prep Apodemia data for analyses
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-21

################################################################################
# SETUP
# Load dependancies & establish data files

# Load dependancies
#install.packages("adegenet")
library("adegenet")
#install.packages("hierfstat")
library("hierfstat")
source(file = "functions/apodemia-functions.R")

# Establish data files
str.file = "data/Apodemia_0.9-noDockweiler-GNPSK.str"
localities.file = "data/Apo_localities.txt"

################################################################################
# GENETIC DATA
# Read in STRUCTURE-formatted file, calculate pairwise Fst matrix, and save the 
# matrix to a file

# Read in file
apo.str.genind <- read.structure(file = str.file,
                                 n.ind = 102,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")

# Calculate pairwise Fst matrix
start <- Sys.time()
pairwise.fst <- genet.dist(apo.str.genind, method = "WC84") # ~10.5 minutes
end <- Sys.time()
total.time <- end - start
cat("Time to calculate pairwise Fst matrix: ", total.time, "\n", sep = "")

# Save file
save(pairwise.fst, file = "output/pairwise-fst.RData")

################################################################################
# GEOGRAPHICAL DATA
# Read in localities file, reconcile it with populations included in STRUCTURE 
# file, calculate geographic distance matrix, and save localities and 
# geographic distance matrix to file



# Read in localities file

# Reconcile localities with populations included in STRUCTURE file

# Calculate geographic distance matrix

# Save localities and geographic distance matrix to file


