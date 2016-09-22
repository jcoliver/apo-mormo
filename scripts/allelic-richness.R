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

genind.file = "output/genind-object.RData"

################################################################################
# DATA PREP
load(file = genind.file) # apo.str.genind

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
