# Apodemia mormo isolation by distance test
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-06

rm(list = ls())
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
# 1. IBD s.s
# 2. Partial Mantel, testing for North-South differentiation beyond IBD
# 3. Partial Mantel, testing for langei differentiation beyond IBD
# 4. IBD s.s. on just North populations
# 5. IBD s.s. on just South populations

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

########################################
# IBD on SOUTHern populations only
south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
south.pop.numbers <- localities$pop.number[which(localities$pop.name %in% south.names)]

# Need to extract part of distance matrices
south.p.fst <- as.matrix(p.fst)[as.character(south.pop.numbers), as.character(south.pop.numbers)]
south.p.fst <- as.dist(south.p.fst)

south.geo.dist <- as.matrix(geo.dist)[as.character(south.pop.numbers), as.character(south.pop.numbers)]
south.geo.dist <- as.dist(south.geo.dist)

ibd <- mantel(xdis = south.geo.dist, 
              ydis = south.p.fst, 
              method = "pearson", 
              permutations = 1000, 
              parallel = 1)
# IBD result
sink(file = ibd.results.file, append = TRUE)
cat("\n--------------------------------------------------------------------------------")
cat("\nIBD South only", sep = "")
print(ibd)
sink()

########################################
# IBD on NORTHern populations only
south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
north.pop.numbers <- localities$pop.number[which(!(localities$pop.name %in% south.names))]

# Need to extract part of distance matrices
north.p.fst <- as.matrix(p.fst)[as.character(north.pop.numbers), as.character(north.pop.numbers)]
north.p.fst <- as.dist(north.p.fst)

north.geo.dist <- as.matrix(geo.dist)[as.character(north.pop.numbers), as.character(north.pop.numbers)]
north.geo.dist <- as.dist(north.geo.dist)

ibd <- mantel(xdis = north.geo.dist, 
              ydis = north.p.fst, 
              method = "pearson", 
              permutations = 1000, 
              parallel = 1)
# IBD result
sink(file = ibd.results.file, append = TRUE)
cat("\n--------------------------------------------------------------------------------")
cat("\nIBD North only", sep = "")
print(ibd)
sink()

########################################
# IBD on NORTHern populations only, partial Mantel on langei
# Partial Mantel to see if langei difference is beyond IBD
south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
north.pop.numbers <- localities$pop.number[which(!(localities$pop.name %in% south.names))]

# Need to extract part of distance matrices
north.p.fst <- as.matrix(p.fst)[as.character(north.pop.numbers), as.character(north.pop.numbers)]
north.p.fst <- as.dist(north.p.fst)

north.geo.dist <- as.matrix(geo.dist)[as.character(north.pop.numbers), as.character(north.pop.numbers)]
north.geo.dist <- as.dist(north.geo.dist)

langei <- "langei"
north.localities <- localities[!c(localities$pop.name %in% south.names), ]
langei.indicator <- IndicatorMatrix(pop.set = langei, localities = north.localities)

ibd <- mantel.partial(xdis = langei.indicator,
              ydis = north.p.fst, 
              zdis = north.geo.dist, 
              method = "pearson", 
              permutations = 1000, 
              parallel = 1)

# IBD result
sink(file = ibd.results.file, append = TRUE)
cat("\n--------------------------------------------------------------------------------")
cat("\nIBD North only, partial Mantel for langei", sep = "")
print(ibd)
sink()
