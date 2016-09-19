# Script for Monmonier analysis of Apodemia data
# Jeff Oliver
# 2016-01-25
# jcoliver@email.arizona.edu

# Much comes from adegenet-tutorial.pdf, which can be downloaded at:
# https://github.com/thibautjombart/adegenet/raw/master/tutorials/tutorial-basics.pdf

setwd("~/Documents/Other Work/Apodemia/Apo-R")

# install.packages("adegenet")
library("adegenet")
# library(hierfstat)
# Read in data
# apo.df <- read.delim("data/filtered_Apo.str", header = FALSE, na.strings = "-9", stringsAsFactors = FALSE)
apo.df <- read.delim("data/filtered_noGNPSK.str", header = FALSE, na.strings = "-9", stringsAsFactors = FALSE)
apo.locs <- read.delim("data/Apo_localities.txt", header = TRUE)

# Drop GrasslandsNPSK
apo.locs <- apo.locs[apo.locs$pop.name != "GrasslandsNPSK",]
apo.locs$pop.name <- factor(apo.locs$pop.name)

# Need to create a vector of population factors...The value in apo.locs$pop.number should 
# match up to value in the second column of apo.df
pop.vector <- factor(x = rep(x = NA, times = length(apo.df[,2])), levels = apo.locs$pop.name)
for (i in 1:length(pop.vector)) {
  pop.number <- apo.df[i, 2]
  pop.name <- apo.locs$pop.name[which(apo.locs$pop.number == pop.number)]
  pop.vector[i] <- pop.name
}
rm(i, pop.name, pop.number)

# Create genid object from data
# apo.genid <- df2genind(X = apo.df[,3:433], ploidy = 2, ind.names = apo.df[,1], pop = apo.df[,2], sep = "")
apo.genid <- df2genind(X = apo.df[,3:433], ploidy = 2, ind.names = apo.df[,1], pop = pop.vector, sep = "")
rm(apo.df, pop.vector)

# We need genetic distance and geographic distance
# Start with Euclidian distance for genetic distance (see p. 69 in adegenet-tutorial.pdf)
genetic.dist <- dist(apo.genid@tab)

# For geographic distance, we need to convert lat/long to UTM and provide:
#   "a list containing: xy"
#   This is really a matrix of doubles
#   Need to convert lat/long to UTM see adegenet-manual.pdf, 
#   under the spca description, see description of xy Argument (p. 54)
#   The convUL function of PBSmapping should be able to do this; see
#   PBSmapping manual p. 29 and example at https://gist.github.com/bmpvieira/3374580

# install.packages("PBSmapping", dependencies = TRUE)
library("PBSmapping")
coord <- data.frame(pop = apo.locs$pop.name, X = apo.locs$longitude, Y = apo.locs$latitude)
attr(coord, "projection") <- "LL"
attr(coord, "zone") <- NA
apo.utm <- convUL(coord)
# apo.utm.mat <- as.matrix(x = apo.utm)
# The xy object in apo.genid needs to have a row for 
# every *individual*
apo.ind.mat <- matrix(data = NA, nrow = length(apo.genid@pop), ncol = 2)
for (i in 1:length(apo.genid@pop)) {
  pop.row <- which(apo.utm$pop == apo.genid@pop[i])
  # It appears (looking at the adgenet example sim2pop data) that each
  # individual requires a unique X,Y coordinate
  apo.ind.mat[i, 1] <- apo.utm$X[pop.row]
  apo.ind.mat[i, 2] <- apo.utm$Y[pop.row]
}
colnames(apo.ind.mat) <- c("x", "y")
# It appears (looking at the adgenet example sim2pop data) that each
# individual requires a unique X,Y coordinate, so jitter points 
# (per adgenet response to chooseCN)
apo.ind.mat[, 1] <- jitter(x = as.vector(apo.ind.mat[, 1]))
apo.ind.mat[, 2] <- jitter(x = as.vector(apo.ind.mat[, 2]))
apo.genid@other$xy <- apo.ind.mat
rm(apo.utm, apo.ind.mat, coord)

# Connectivity network
gab <- chooseCN(xy = apo.genid@other$xy, ask = FALSE, type = 2)

# Get rid of random noise
pco1 <- dudi.pco(d = genetic.dist, scannf = FALSE, nf = 1)
barplot(pco1$eig, main = "Eigenvalues of Genetic Distance PCO")
genetic.dist.eig <- dist(pco1$li)

mon1 <- monmonier(xy = apo.genid@other$xy, dist = genetic.dist.eig, cn = gab)
# Chose default threshold

# Need to do better colors
plot(mon1, add.arrows = FALSE, bwd = 10, col = "black")
points(apo.genid@other$xy, cex = 2, pch = 20,
       col = fac2col(pop(apo.genid), col.pal = spectral))
legend("topright", legend = levels(pop(apo.genid)), pch = c(20), 
       col = spectral(length(levels(pop(apo.genid)))), pt.cex = 2)


# Trying two runs
mon2 <- monmonier(xy = apo.genid@other$xy, dist = genetic.dist.eig, cn = gab, nrun = 2)
plot(mon2, add.arrows = FALSE, bwd = 10, col = c("red", "magenta", "blue", "black", "cyan", "green"))

