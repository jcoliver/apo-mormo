# Rarefaction test 3: Work on genotype format data
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-03

rm(list = ls())

source(file = "functions/rarefaction.R")

################################################################################
# Can create N.matrix from using table on a data frame
locus.data <- data.frame(pop = c(rep(x = 1, times = 6), rep(x = 2, times = 6), rep(x = 3, times = 4), rep(x = 4, times = 4)),
                         allele = c(rep(x = "10", times = 3), rep(x = "11", times = 3), rep(x = "10", times = 6), rep(x = "10", times = 4), rep(x = "12", times = 4)))

# Transpose to get N matrix
N.matrix <- t(unclass(table(locus.data)))

################################################################################
rich.g2 <- calcRichness.all(N = N.matrix, g = 2)
rich.g4 <- calcRichness.all(N = N.matrix, g = 4)
private.g2 <- calcPrivate.all(N = N.matrix, g = 2)
private.g4 <- calcPrivate.all(N = N.matrix, g = 4)

################################################################################
pass <- all(rich.g2 == c(1.6, 1, 1, 1)) & 
  all(rich.g4 == c(2, 1, 1, 1)) &
  all(private.g2 == c(0.8, 0, 0, 1)) &
  all(private.g4 == c(1, 0, 0, 1))

cat("Test 3 pass = ", pass, "\n", sep = "")

rm(list = ls())
