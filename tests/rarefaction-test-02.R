# Rarefaction test 2: three alleles, four populations
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-03

rm(list = ls())

source(file = "functions/rarefaction.R")

################################################################################
# Set up data to use
N.matrix <- matrix(data = c(3, 6, 4, 0, 
                            3, 0, 0, 0, 
                            0, 0, 0, 4), nrow = 3, ncol = 4, byrow = TRUE)

################################################################################
private.g2 <- calcPrivate.all(N = N.matrix, g = 2)
private.g4 <- calcPrivate.all(N = N.matrix, g = 4)

################################################################################
pass <- all(private.g2 == c(0.8, 0, 0, 1)) &
  all(private.g4 == c(1, 0, 0, 1))

cat("Test 2 pass = ", pass, "\n", sep = "")

rm(list = ls())
