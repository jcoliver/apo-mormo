# Rarefaction test 9: Move loci-loop to function
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-10

rm(list = ls())

library("adegenet")

source(file = "functions/rarefaction-functions.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

################################################################################
# Identify those populations were are interested in
pops <- c(1, 9, 11, 12, 13)
pop.rows <- which(as.integer(as.character(apo.str.genind@pop)) %in% pops)

# Subset data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[pop.rows, 1:20]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:20])
test.data$pop <- factor(apo.str.genind@pop[pop.rows])

rarefied.values <- rarefiedMatrices(data = test.data, g = 6)
richness.matrix <- rarefied.values$richness
private.matrix <- rarefied.values$private

################################################################################
richness.test <- c(1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 
                   NA, NA, NA, NA, 2, 1, 2, 1, 1, 1,
                   1, 1, 1, 1, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
################################################################################
richness.matrix <- round(x = richness.matrix, digits = 3)
richness.test <- round(x = richness.test, digits = 3)

pass <- all(as.vector(richness.matrix) == richness.test, na.rm = TRUE) &
  all(is.na(private.matrix))

cat("Test 9 pass = ", pass, "\n", sep = "")

rm(list = ls())