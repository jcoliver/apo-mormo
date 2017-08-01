# Prep Apodemia data for analyses
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-21

################################################################################
# SUMMARY
# * Reads in genetic data from a STRUCTURE-formatted file and calculates pairwise
# Fst matrix. Saves the matrix to output folder
# * Reads in locality data file, which includes latitude / longitude coordinates,
# and reconciles it with the populations included in the STRUCTURE file. Then 
# calculates pairwise geographic distance matrix and saves reconciled localities
# and geographic distance matrix to output folder

################################################################################
# SETUP
# Load dependancies & establish data and output files

# Load dependancies
#install.packages("adegenet")
library("adegenet")
#install.packages("hierfstat")
library("hierfstat")
source(file = "functions/apodemia-functions.R")

# Establish data files
str.file = "data/Apodemia_0.9-noGNPSK.str"
num.inds <- 107 # 102
genind.out.file = "output/genind-object.RData"
fst.out.file = "output/pairwise-fst.RData"
localities.file = "data/Apo_localities.txt"
localities.out.file = "output/reconciled-localities.RData"
geo.dist.out.file = "output/pairwise-geo-dist.RData"
omit.locs <- c("GrasslandsNPSK") # c("Dockweiler", "GrasslandsNPSK")

################################################################################
# GENETIC DATA
# Read in STRUCTURE-formatted file, calculate pairwise Fst matrix, and save the 
# matrix to a file

# Read in file
apo.str.genind <- read.structure(file = str.file,
                                 n.ind = num.inds,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")

save(apo.str.genind, file = genind.out.file)

# Calculate pairwise Fst matrix
start <- Sys.time()
pairwise.fst <- genet.dist(apo.str.genind, method = "WC84") # ~10.5 minutes
end <- Sys.time()
total.time <- end - start
cat("Time to calculate pairwise Fst matrix: ", total.time, "\n", sep = "")

# Save file
save(pairwise.fst, file = fst.out.file)

################################################################################
# GEOGRAPHICAL DATA
# Reconcile localities in file with populations included in STRUCTURE file, 
# calculate geographic distance matrix, and save localities and geographic 
# distance matrix to file

# Reconcile localities with populations included in STRUCTURE file
localities <- FormatLocalities(file = localities.file,
                               genind = apo.str.genind, 
                               omit = omit.locs)

# Calculate geographic distance matrix
geo.dist <- GeoDistances(localities = localities)

# Save localities and geographic distance matrix to file
save(localities, file = localities.out.file)
save(geo.dist, file = geo.dist.out.file)

