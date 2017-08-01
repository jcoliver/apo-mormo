# Rarefaction test 10: Work on rarefaction for pops with too low sampled genes
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-10

rm(list = ls())

library("adegenet")

source(file = "functions/rarefaction-functions.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

################################################################################
# DESCRIPTION OF ISSUE:
# Example:
# with populations 1, 9, 11, 12, 13 and loci 1-10, g = 6
# population 13 has 0 samples for loci 5-10
# Richness has NA for population 13 loci 5-10
# Richness has NA for ALL POPULATIONS in loci 5-10
#
# What causes NAs in the output matrices?
# g is greater than # of samples
# In richness, manifests as NA for any population with fewer sampled genes than 
# g; this seems like it is OK behavior
# In private, manifests as NA for *all* populations for a locus with fewer 
# sampled genes in *any* population. This is NOT OK. Need to remove the 
# offending populations before private allele calculations proceed, but will 
# have to keep track of which have been removed so the resultant matrix is 
# still indexed correctly

################################################################################
# Identify those populations were are interested in
pops <- c(1, 9, 11, 12, 13)
pop.rows <- which(as.integer(as.character(apo.str.genind@pop)) %in% pops)

# Subset data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[pop.rows, 1:20]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:20])
test.data$pop <- factor(apo.str.genind@pop[pop.rows])

rarefied.values <- rarefiedMatrices(data = test.data, g = 14)
richness.matrix <- rarefied.values$richness
private.matrix <- rarefied.values$private

################################################################################
richness.test <- c(1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 
                   NA, NA, NA, NA, 2, 1, 2, 1, 1, 1,
                   1, 1, 1, 1, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
private.test <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 
                  NA, NA, NA, NA, 1, 0, 1, 0, 0, 1,
                  0, 0, 0, 0, NA, NA, NA, NA, NA, NA,
                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
################################################################################
richness.matrix <- round(x = richness.matrix, digits = 3)
richness.test <- round(x = richness.test, digits = 3)
private.matrix <- round(x = private.matrix, digits = 3)
private.test <- round(x = private.test, digits = 3)

pass <- all(as.vector(richness.matrix) == richness.test, na.rm = TRUE) &
  all(as.vector(private.matrix) == private.test, na.rm = TRUE)

cat("Test 10 pass = ", pass, "\n", sep = "")

rm(list = ls())