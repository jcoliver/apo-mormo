# Rarefaction test 4: Real data, extracted from genind object (2 pops, 1 locus)
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-03

rm(list = ls())

library("dplyr")
library("adegenet")

source(file = "functions/rarefaction-functions.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

################################################################################
# Grab some data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[1:15, 1:20]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:20])
test.data$pop <- factor(apo.str.genind@pop[1:15])

# Get allele x pop matrix for locus 9
locus.data <- as.data.frame(test.data$tab[, 17:18])
locus.data$pop <- test.data$pop
locus.fac <- factor(test.data$loc.fac[17:18])
locus.id <- unique(test.data$loc.fac[17:18])

# Get counts for each column (allele) of the data
allele.counts <- locus.data %>% 
  group_by(pop) %>% 
  summarize_all(funs(sum), na.rm = TRUE)

# Pull out counts (first column is pop id)
count.matrix <- as.matrix(allele.counts[, c(2:ncol(allele.counts))])

# Set row names from values of pop (odd syntax to extract one column 
# from the allele.counts tibble)
rownames(count.matrix) <- as.character(allele.counts[[1]])
colnames(count.matrix) <- gsub(pattern = paste0(as.character(locus.id), "."),
                               replacement = "", 
                               x = colnames(count.matrix))

N.matrix <- t(count.matrix)

################################################################################
richness <- calcRichness.all(N = N.matrix, g = 14)
private <- calcPrivate.all(N = N.matrix, g = 14)

################################################################################
pass <- all(richness == c(2, 2)) & 
  all(private == c(0, 0))

cat("Test 4 pass = ", pass, "\n", sep = "")

rm(list = ls())
