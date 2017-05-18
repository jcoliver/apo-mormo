# Rarefaction testing
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-05-17

################################################################################
rm(list = ls())

source(file = "functions/rarefaction.R")

## TESTS
N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)

calcRichness(N = N.matrix, j = 1, g = 2)
calcRichness(N = N.matrix, j = 2, g = 2)
calcRichness(N = N.matrix, j = 1, g = 4)
calcRichness(N = N.matrix, j = 2, g = 4)

# P112 = 1 - Q112
calcP(Nj = sum(N.matrix[, 1]), Nij = N.matrix[1, 1], g = 2)

# Q122
calcQ(Nj = sum(N.matrix[, 2]), Nij = N.matrix[1, 2], g = 2)

# P212 = 1 - Q212
calcP(Nj = sum(N.matrix[, 1]), Nij = N.matrix[2, 1], g = 2)

# Q222
calcQ(Nj = sum(N.matrix[, 2]), Nij = N.matrix[2, 2], g = 2)

# Prob allele *not* found
calcQ.v(N.col = N.matrix[, 2], g = 2)
calcQ.v(N.col = N.matrix[, 1], g = 2)

# Prob allele found
calcP.v(N.col = N.matrix[, 2], g = 2)
calcP.v(N.col = N.matrix[, 1], g = 2)

# Prob for single allele **doesn't work if N.col = 0 (for allele not found in population)
# CANNOT work, because Nj (total # of genes surveyed from population is not included)
calcP.v(N.col = N.matrix[1, 2], g = 2)
calcP.v(N.col = N.matrix[2, 2], g = 2)

# Allele richness for individual populations
calcRichness.v(N.col = N.matrix[, 1], 4)
calcRichness.v(N.col = N.matrix[, 2], 4)

# Allele richness for each population
calcRichness.all(N = N.matrix, g = 2)
calcRichness.all(N = N.matrix, g = 2, loop = FALSE)
calcRichness.all(N = N.matrix, g = 4)
calcRichness.all(N = N.matrix, g = 4, loop = FALSE)
