# Rarefaction test 12: Real data, ~ 1000 loci
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-11

rm(list = ls())

library("adegenet")

source(file = "functions/rarefaction.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

################################################################################
# Identify those populations were are interested in
pops <- c(1, 9, 11, 12, 13)
pop.rows <- which(as.integer(as.character(apo.str.genind@pop)) %in% pops)

# Subset data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[pop.rows, 1:1925]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:1925])
test.data$pop <- factor(apo.str.genind@pop[pop.rows])

rarefied.values <- rarefiedMatrices(data = test.data, 
                                    g = 10, 
                                    display.progress = TRUE)
richness.matrix <- rarefied.values$richness
private.matrix <- rarefied.values$private

################################################################################
richness.test <- 4240.792
private.test <- 396.6959
################################################################################
richness.matrix.sum <- round(x = sum(x = richness.matrix, na.rm = TRUE), 
                             digits = 3)
richness.test <- round(x = richness.test, digits = 3)
private.matrix.sum <- round(x = sum(x = private.matrix, na.rm = TRUE), 
                            digits = 3)
private.test <- round(x = private.test, digits = 3)

pass <- richness.matrix.sum == richness.test &&
  private.matrix.sum == private.test

cat("Test 12 pass = ", pass, "\n", sep = "")

rm(list = ls())
