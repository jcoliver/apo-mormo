# Rarefaction test 1: two alleles, two populations
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-03

rm(list = ls())

source(file = "functions/rarefaction-functions.R")

################################################################################
# Set up data to use
N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)

################################################################################
# Allele richness for each population
rich.g2 <- calcRichness.all(N = N.matrix, g = 2)
rich.g4 <- calcRichness.all(N = N.matrix, g = 4)

# Private alleles for all populations
private.g2 <- calcPrivate.all(N = N.matrix, g = 2)
private.g4 <- calcPrivate.all(N = N.matrix, g = 4)

################################################################################
pass <- all(rich.g2 == c(1.5, 1)) & 
  all(rich.g4 == c(2, 1)) &
  all(private.g2 == c(0.5, 0)) &
  all(private.g4 == c(1, 0))

cat("Test 1 pass = ", pass, "\n", sep = "")

rm(list = ls())