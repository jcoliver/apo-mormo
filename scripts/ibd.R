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

########################################
# IBD
ibd <- mantel(xdis = geo.dist, 
              ydis = p.fst, 
              method = "pearson", 
              permutations = 1000, 
              parallel = 1)

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



localities <- FormatLocalities(file = "data/Apo_localities.txt",
                               genind = apo.str.genind, 
                               omit = c("Dockweiler", "GrasslandsNPSK"))

geo.dist <- GeoDistances(localities = localities)

# Mantel tests
# Want Fst/(1 - Fst) as our differentiation matrix
p.fst <- pairwise.fst/(1 - pairwise.fst)
# Want log-transformed geographic distances, per Rousset 1997
# Rousset does not explicitly mention the base of the log function, 
# but http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-6-13
# say base-10, so we go with that
geo.dist <- log(x = geo.dist, base = 10)

ibd <- mantel(xdis = geo.dist, ydis = p.fst, method = "pearson", permutations = 1000, parallel = 1) # parallel = getOption("mc.cores")
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = geo.dist, ydis = p.fst, method = "pearson", permutations = 1000,      parallel = 1) 
# 
# Mantel statistic r: 0.4995 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.191 0.256 0.303 0.366 
# Permutation: free
# Number of permutations: 1000

# Partial Mantel to see if North-South difference is beyond IBD
south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
north.south <- IndicatorMatrix(pop.set = south.names, localities = localities)

ibd.ns <- mantel.partial(xdis = north.south, ydis = p.fst, zdis = geo.dist, method = "pearson", permutations = 1000, parallel = 1)
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = north.south, ydis = p.fst, zdis = geo.dist, method = "pearson", permutations = 1000, parallel = 1) 
# 
# Mantel statistic r: -0.1584 
#       Significance: 0.82018 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.200 0.271 0.324 0.358 
# Permutation: free
# Number of permutations: 1000

# Partial Mantel to see if langei difference is beyond IBD
langei <- "langei"
langei.indicator <- IndicatorMatrix(pop.set = langei, localities = localities)
ibd.langei <- mantel.partial(xdis = langei.indicator, ydis = p.fst, zdis = geo.dist, method = "pearson", permutations = 1000, parallel = 1)
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = langei.indicator, ydis = p.fst, zdis = geo.dist,      method = "pearson", permutations = 1000, parallel = 1) 
# 
# Mantel statistic r: 0.4635 
#       Significance: 0.15584 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.464 0.552 0.552 0.552 
# Permutation: free
# Number of permutations: 1000
