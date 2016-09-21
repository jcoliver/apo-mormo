# Apodemia mormo isolation by distance test
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-06

################################################################################
# SUMMARY
# * Reads in fst and geographic distance matrices, as well as localities object
# * Performs isolation by distance Mantel test
# * Performs partial Mantel test for differentiation between north and south 
# populations after conrolling for IBD
# * Performs partial Mantel test for differentiation of langei from remaining 
# populations after controlling for IBD

################################################################################
# SETUP
# Load dependancies & establish data files

# Load dependancies
#install.packages("vegan")
library("vegan")
source(file = "functions/apodemia-functions.R")

# Establish data files
pairwise.fst.file = "output/pairwise-fst.RData"
localities.file = "output/reconciled-localities.RData"
geo.dist.file = "output/pairwise-geo-dist.RData"
ibd.results.file = "output/ibd-results.txt"

################################################################################
# DATA PREP

# Load pairwise Fst matrix
load(file = pairwise.fst.file)
# Want Fst/(1 - Fst) as our differentiation matrix
p.fst <- pairwise.fst/(1 - pairwise.fst)

# Load localities
load(file = localities.file)

# Load geographic distance matrix
load(file = geo.dist.file)
# Want log-transformed geographic distances, per Rousset 1997
# Rousset does not explicitly mention the base of the log function, 
# but http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-6-13
# say base-10, so we go with that
geo.dist <- log(x = geo.dist, base = 10)

################################################################################
# IBD ANALYSES

# Setup results file
sink(file = ibd.results.file, append = FALSE)
cat("IBD Results", "\n", sep = "")
cat(format.Date(Sys.Date(), "%Y-%m-%d"), "\n", sep = "")
cat("Fst file: ", pairwise.fst.file, "\n", sep = "")
cat("Localities file: ", pairwise.fst.file, "\n", sep = "")
cat("Geographic distances file: ", geo.dist.file, "\n", sep = "")
sink()

########################################
# IBD
ibd <- mantel(xdis = geo.dist, 
              ydis = p.fst, 
              method = "pearson", 
              permutations = 1000, 
              parallel = 1)
# IBD result
sink(file = ibd.results.file, append = TRUE)
cat("\n--------------------------------------------------------------------------------")
cat("\nIBD Mantel", sep = "")
print(ibd)
sink()

########################################
# Partial Mantel to see if North-South difference is beyond IBD
south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
north.south <- IndicatorMatrix(pop.set = south.names, localities = localities)

ibd.ns <- mantel.partial(xdis = north.south, 
                         ydis = p.fst, 
                         zdis = geo.dist, 
                         method = "pearson", 
                         permutations = 1000, 
                         parallel = 1)

sink(file = ibd.results.file, append = TRUE)
cat("\n--------------------------------------------------------------------------------")
cat("\nIBD partial Mantel North-South", sep = "")
print(ibd.ns)
sink()

########################################
# Partial Mantel to see if langei difference is beyond IBD
langei <- "langei"
langei.indicator <- IndicatorMatrix(pop.set = langei, localities = localities)
ibd.langei <- mantel.partial(xdis = langei.indicator, 
                             ydis = p.fst, 
                             zdis = geo.dist, 
                             method = "pearson", 
                             permutations = 1000, 
                             parallel = 1)

sink(file = ibd.results.file, append = TRUE)
cat("\n--------------------------------------------------------------------------------")
cat("\nIBD partial Mantel langei", sep = "")
print(ibd.langei)
sink()
