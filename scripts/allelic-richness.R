# Script for Apodemia allelic richness analyses
# Jeff Oliver
# 2016-08-23
# jcoliver@email.arizona.edu

# INCOMPLETE

################################################################################
# SUMMARY


################################################################################
# SETUP
#install.packages("hierfstat")
library("hierfstat")
#install.packages("tidyr")
library("tidyr")

source(file = "functions/rarefaction.R")

genind.file = "output/genind-object.RData"

min.pop.size <- 4

################################################################################
# DATA PREP
load(file = genind.file) # apo.str.genind

# Need to get data in proper format...
richness <- matrix(data = 0, nrow = length(levels(apo.str.genind@loc.fac)),
                   ncol = length(levels(apo.str.genind@pop)))
rownames(richness) <- levels(apo.str.genind@loc.fac)
colnames(richness) <- levels(apo.str.genind@pop)

# Looping for now, but eventually use an apply function
# Loop over each locus
for (locus in levels(apo.str.genind@loc.fac)) {
  # prep data for calcPrivate.all and calcRichness.all
  N.matrix
  # do calculations for all pops for that one locus
  locus.richness <- calcRichness.all(N = N.matrix, g = min.pop.size)
  # store results for that locus in richness matrix
  richness[rownames(richness) == locus, ] <- locus.richness
}

# Old, error-laden way:

# Data frame for allelic richness
apo.df <- cbind(pop = apo.str.genind@pop, apo.str.genind@tab)
rownames(apo.df) <- rownames(apo.str.genind@tab)
richness <- allelic.richness(data = apo.df, diploid = TRUE)



# Old (original) way
apo.df <- read.delim(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     na.strings = "-9", 
                     header = FALSE)
# Need to find population with fewest individuals for rarefaction
min.count <- min(table(apo.df[, 2]))
min.count <- min.count * 2 # for diploids

richness <- allelic.richness(data = apo.df[, -1], min.n = min.count, diploid = TRUE)

richness.df <- data.frame(richness$Ar)
colnames(richness.df) <- paste0("pop.", unique(apo.df[, 2]))
mean.richness <- colMeans(x = richness.df, na.rm = TRUE)

richness.long <- gather(data = richness.df, key = population, value = richness)
