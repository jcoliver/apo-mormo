# Script for Apodemia allelic richness analyses
# Jeff Oliver
# 2016-08-23
# jcoliver@email.arizona.edu

#install.packages("hierfstat")

apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", header = FALSE, na.strings = "-9", stringsAsFactors = FALSE)

# apo.df is:
# First column: sample name
# Second column: population (numeric)
# Remaining column: SNP data, with two columns per locus

# Need to find population with fewest individuals for rarefaction
min.count <- min(table(apo.df[, 2]))
min.count <- min.count * 2 # for diploids

library("hierfstat")
richness <- allelic.richness(data = apo.df[, -1], min.n = min.count, diploid = TRUE)

richness.df <- data.frame(richness$Ar)
colnames(richness.df) <- paste0("pop.", unique(apo.df[, 2]))
mean.richness <- colMeans(x = richness.df, na.rm = TRUE)

#install.packages("tidyr")
library("tidyr")
richness.long <- gather(data = richness.df, key = population, value = richness)

# South populations: PointLoma, WildhorseMeadows, CampPendleton, Borrego
