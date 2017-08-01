# Rarefaction test 13: Real data, All populations, ~ 1000 loci
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-12

rm(list = ls())

library("adegenet")

source(file = "functions/rarefaction-functions.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

################################################################################
# Subset data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[, 1:1925]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:1925])
test.data$pop <- factor(apo.str.genind@pop)

rarefied.values <- rarefiedMatrices(data = test.data, 
                                    g = 10, 
                                    display.progress = TRUE)
richness.matrix <- rarefied.values$richness
private.matrix <- rarefied.values$private

################################################################################
richness.test <- 12718.92
private.test <- 393.4783
################################################################################
richness.matrix.sum <- round(x = sum(x = richness.matrix, na.rm = TRUE), 
                             digits = 2)
richness.test <- round(x = richness.test, digits = 2)
private.matrix.sum <- round(x = sum(x = private.matrix, na.rm = TRUE), 
                            digits = 2)
private.test <- round(x = private.test, digits = 2)

pass <- richness.matrix.sum == richness.test &&
  private.matrix.sum == private.test

cat("Test 13 pass = ", pass, "\n", sep = "")

rm(list = ls())
