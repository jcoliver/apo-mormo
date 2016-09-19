# Script for Apodemia IBD analyses
# Jeff Oliver
# 2016-01-25
# jcoliver@email.arizona.edu

#install.packages("hierfstat")
#install.packages("adegenet")

# OR, just load(file = "output/Apo_pairwise_Fst.RData")
# apo.df <- read.delim("data/filtered_Apo.str", header = FALSE, na.strings = "-9", stringsAsFactors = FALSE)
apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", header = FALSE, na.strings = "-9", stringsAsFactors = FALSE)

# apo.df is:
# First column: sample name
# Second column: population (numeric)
# Remaining column: SNP data, with two columns per locus

# See https://nescent.github.io/popgenInfo/DifferentiationSNP.html
library(adegenet)
library(hierfstat)
apo.genid <- df2genind(X = apo.df[,3:433], ploidy = 2, ind.names = apo.df[,1], pop = apo.df[,2], sep = "")
# Calculate pairwise Fst (takes ~3 minutes)
p.fst <- genet.dist(apo.genid, method="WC84") # A 12 x 12 matrix (lower triangle only)
# Much more reasonable values... mean: 0.03 (0.001 - 0.06)
save(p.fst, file = "output/Apo_pairwise_Fst.RData")

load(file = "output/Apo_pairwise_Fst.RData")

# For isolation by distance, need a distance matrix among populations.
# Need a data.frame with lat/long (decimal degrees) and population numbers
# pop.name  pop.number  latitude  longitude
localities <- read.delim(file = "data/Apo_localities.txt", header = TRUE)

# Take out Dockweiler & Grasslands populations
localities <- localities[-which(localities$pop.name %in% c("Dockweiler", "GrasslandsNPSK")), ]

# Need to re-number populations in localities to match their numbers in apo.df
apo.localities <- apo.df[, c(1, 2)] # First column of data frame has sample name (which includes pop)
colnames(apo.localities) <- c("pop.name", "pop.number")
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

# Convert from decimal degrees to radians
localities$lat.rad <- (localities$latitude * pi)/180
localities$long.rad <- (localities$longitude * pi)/180

# Use the Haversine formula for great circle distances
# see http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
haversine <- function(lat1, long1, lat2, long2) {
  earth.rad <- 6371 # km
  delta.lat <- lat1 - lat2
  delta.long <- long1 - long2
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = earth.rad * c
  return(d) # Distance in km
}

#geo.dist <- matrix(data = NA, nrow = nrow(localities) - 1, ncol = nrow(localities) - 1)
#colnames(geo.dist) <- localities$pop.number[1:length(localities$pop.number) - 1]
#rownames(geo.dist) <- localities$pop.number[2:length(localities$pop.number)]
geo.dist <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
colnames(geo.dist) <- localities$pop.number[1:length(localities$pop.number)]
rownames(geo.dist) <- localities$pop.number[1:length(localities$pop.number)]

for (i in 1:(length(localities$lat.rad) - 1)) {
  lat1 <- localities$lat.rad[i]
  long1 <- localities$long.rad[i]
  for (j in (i + 1):length(localities$lat.rad)) {
    lat2 <- localities$lat.rad[j]
    long2 <- localities$long.rad[j]
    d <- haversine(lat1 = lat1, long1 = long1, lat2 = lat2, long2 = long2)
    geo.dist[j, i] <- d
  }
}

rm(d, i, j, lat1, long1, lat2, long2)

geo.dist <- as.dist(geo.dist)
# Checking work NOTE these are almost certainly wrong 2016-09-15
# 1 to 4 (14?) 19.94 km
# 8 to 15 727.95 km
#haversine(localities$lat.rad[1], localities$long.rad[1], localities$lat.rad[4], localities$long.rad[4])
#haversine(localities$lat.rad[8], localities$long.rad[8], localities$lat.rad[15], localities$long.rad[15])


# Mantel tests
# Want Fst/(1 - Fst) as our differentiation matrix
p.fst <- p.fst/(1 - p.fst)
# Want log-transformed geographic distances, per Rousset 1997
# Rousset does not explicitly mention the base of the log function, 
# but http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-6-13
# say base-10, so we go with that
geo.dist <- log(x = geo.dist, base = 10)

# install.packages("vegan")
library(vegan)
ibd <- mantel(xdis = geo.dist, ydis = p.fst, method = "pearson", permutations = 1000, parallel = 1) # parallel = getOption("mc.cores")
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = geo.dist, ydis = p.fst, method = "pearson", permutations = 1000,      parallel = 1) 
# 
# Mantel statistic r:  0.79 
#       Significance: 0.000999 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.182 0.242 0.306 0.363 
# Permutation: free
# Number of permutations: 1000

# Try testing for additional differentiation in langei
lang.diff <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
lang.row <- which(localities$pop.name == "langei")
for (i in 1:(length(localities$lat.rad) - 1)) {
  for (j in (i + 1):length(localities$lat.rad)) {
    d <- 0
    if (i == lang.row || j == lang.row) {
      d <- 1
    }
    lang.diff[j, i] <- d
  }
}
rm(i, j, d)
lang.diff <- as.dist(lang.diff)

# Now see if langei are different even when controlling for isolation by distance
ibd.par <- mantel.partial(xdis = lang.diff, ydis = p.fst, zdis = geo.dist, method = "pearson", permutations = 1000, parallel = 1)
# Partial Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel.partial(xdis = lang.diff, ydis = p.fst, zdis = geo.dist,      method = "pearson", permutations = 1000, parallel = 1) 
# 
# Mantel statistic r:  0.38 
#       Significance: 0.082917 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.358 0.380 0.380 0.380 
# Permutation: free
# Number of permutations: 1000

# Plotting distance vs. Fst (both transformed)
# Coloring comparisons with langei red
cols <- lang.diff
cols[cols == 0] <- "black"
cols[cols == 1] <- "red"
plot(x = geo.dist, y = p.fst, col = cols, pch = 19, xlab = "Log(distance)", ylab = "Fst/1 - Fst")
abline(lm(p.fst ~ geo.dist))


# Partial Mantel to see if North-South difference is beyond IBD
#south.names <- c("PointLoma", "WildhorseMeadows", "Dockweiler", "CampPendleton", "GrasslandsNPSK", "Borrego")
south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
south <- which(localities$pop.name %in% south.names)
north.south <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
for (i in 1:length(localities$lat.rad) - 1) {
  for (j in 1:length(localities$lat.rad)) {
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
# mantel.partial(xdis = north.south, ydis = p.fst, zdis = geo.dist,      method = "pearson", permutations = 1000, parallel = 1) 
# 
# Mantel statistic r: 0.2565 
#       Significance: 0.032967 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.196 0.234 0.263 0.310 
# Permutation: free
# Number of permutations: 1000


# Plotting distance vs. Fst (both transformed)
# Coloring comparisons with North vs. South red
cols <- north.south
cols[cols == 0] <- "white"
cols[cols == 1] <- "black"
date.filename <- format(Sys.Date(), "%Y-%m-%d")
pdf(file = paste0("output/Apo-ibd-graph-", date.filename, ".pdf"), useDingbats = FALSE)
  plot(x = geo.dist, y = p.fst, col = "black", pch = 21, bg = cols, xlab = "Log(distance)", ylab = "Fst/1 - Fst")
# abline(lm(p.fst ~ geo.dist))
dev.off()

# Doing some clustering, we see two groups
# See https://nescent.github.io/popgenInfo/DifferentiationSNP.html
grp <- find.clusters(apo.genid, max.n.clust = 10, n.pca = 20, choose.n.clust = FALSE) 
dapc1 <- dapc(apo.genid, grp$grp, n.pca = 20, n.da = 6) 
scatter(dapc1)

