# Apodemia mormo isolation by distance test
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-06

#install.packages("hierfstat")
#install.packages("adegenet")
library("adegenet")
library("hierfstat")

################################################################################
# Load structure file to genind object
apo.str.genind <- read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str",
                                 n.ind = 102,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")
#' /// GENIND OBJECT /////////
#'   
#'   // 102 individuals; 4,057 loci; 7,808 alleles; size: 4.8 Mb
#' 
#' // Basic content
#' @tab:  102 x 7808 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 1-2)
#' @loc.fac: locus factor for the 7808 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str", 
#'                       n.ind = 102, n.loc = 4057, onerowperind = TRUE, col.lab = 1, 
#'                       col.pop = 2, col.others = 0, row.marknames = 0, NA.char = "-9", 
#'                       sep = "\t")
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 4-12)
#' @other: a list containing: X 

apo.str.summary <- summary(apo.str.genind)
start <- Sys.time()
pairwise.fst <- genet.dist(apo.str.genind, method = "WC84") # ~10.5 minutes
end <- Sys.time()
total.time <- end - start
total.time
save(pairwise.fst, file = "output/pairwise-fst-Apodemia_0.9-noDockweiler-GNPSK.RData")


source(file = "functions/apodemia-functions.R")

localities <- FormatLocalities(file = "data/Apo_localities.txt",
                               genind = apo.str.genind, 
                               omit = c("Dockweiler", "GrasslandsNPSK"))

geo.dist <- GeoDistances(localities = localities)

  


# For isolation by distance, need a distance matrix among populations.
# Need a data.frame with lat/long (decimal degrees) and population numbers
# pop.name  pop.number  latitude  longitude
localities <- read.delim(file = "data/Apo_localities.txt", header = TRUE)

# Take out Dockweiler & Grasslands populations
localities <- localities[-which(localities$pop.name %in% c("Dockweiler", "GrasslandsNPSK")), ]

# Need to re-number populations in localities to match their numbers in apo.df
apo.localities <- data.frame(pop.name = rownames(apo.str.genind@tab), pop.number = apo.str.genind@pop)
apo.localities$pop.name <- gsub("_[0-9A-Za-z]+", "", apo.localities$pop.name)
apo.localities <- apo.localities[order(apo.localities$pop.name), ]
apo.localities <- apo.localities[match(unique(apo.localities$pop.name), apo.localities$pop.name), ]

# Sanity check
if (nrow(localities) != nrow(apo.localities)) {
  stop("Different number of localities in data frames")
}

# Drop the pop.number column, as we'll replace it with a merge
localities <- localities[, -which(colnames(localities) == "pop.number")]
localities <- merge(x = localities, y = apo.localities, by = "pop.name")

# MUST have these in same order as p.fst matrix...could do this better
localities <- localities[order(localities$pop.number), ]
rownames(localities) <- NULL

geo.dist <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
colnames(geo.dist) <- localities$pop.number[1:length(localities$pop.number)]
rownames(geo.dist) <- localities$pop.number[1:length(localities$pop.number)]

for (i in 1:(length(localities$latitude) - 1)) {
  lat1 <- localities$latitude[i]
  long1 <- localities$longitude[i]
  for (j in (i + 1):length(localities$latitude)) {
    lat2 <- localities$latitude[j]
    long2 <- localities$longitude[j]
    d <- Haversine(lat1 = lat1, long1 = long1, lat2 = lat2, long2 = long2)
    geo.dist[j, i] <- d
  }
}

rm(d, i, j, lat1, long1, lat2, long2)

geo.dist <- as.dist(geo.dist)
# Checking work NOTE: indexes do not equal population number!!!
# HullMtn to Borrego
# Haversine(localities$latitude[1], localities$longitude[1], localities$latitude[4], localities$longitude[4])
# 915.8865
# HullMtn to Ladoga
# Haversine(localities$latitude[1], localities$longitude[1], localities$latitude[6], localities$longitude[6])
# 76.89984

# Mantel tests
# Want Fst/(1 - Fst) as our differentiation matrix
p.fst <- pairwise.fst/(1 - pairwise.fst)
# Want log-transformed geographic distances, per Rousset 1997
# Rousset does not explicitly mention the base of the log function, 
# but http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-6-13
# say base-10, so we go with that
geo.dist <- log(x = geo.dist, base = 10)

# install.packages("vegan")
library("vegan")
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
south <- which(localities$pop.name %in% south.names)
north.south <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
for (i in 1:nrow(localities) - 1) {
  for (j in 1:nrow(localities)) {
    d <- 0
    if ((i %in% south && !(j %in% south)) 
        || (j %in% south && !(i %in% south))){
      d <- 1
    }
    north.south[j, i] <- d
  }
}
rm(i, j, d)
north.south <- as.dist(north.south)
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

# Plotting distance vs. Fst (both transformed)
# Coloring comparisons with North vs. South black
cols <- north.south
cols[cols == 0] <- "white"
cols[cols == 1] <- "black"
date.filename <- format(Sys.Date(), "%Y-%m-%d")
pdf(file = paste0("output/Apo-ibd-ns-graph-", date.filename, ".pdf"), useDingbats = FALSE)
plot(x = geo.dist, y = p.fst, col = "black", pch = 21, bg = cols, xlab = "Log(distance)", ylab = "Fst/1 - Fst")
abline(lm(p.fst ~ geo.dist))
dev.off()

