# Rarefaction test 6: Real data, extracted from genind object (4 pops, 10 loci)
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
# Subset data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[1:35, 1:20]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:20])
test.data$pop <- factor(apo.str.genind@pop[1:35])

# Establish matrices that will hold the final results
richness.matrix <- matrix(data = 0, 
                          nrow = length(levels(test.data$loc.fac)),
                          ncol = length(levels(test.data$pop)))

private.matrix <- matrix(data = 0, 
                         nrow = length(levels(test.data$loc.fac)),
                         ncol = length(levels(test.data$pop)))

colnames(richness.matrix) <- colnames(private.matrix) <- as.character(levels(test.data$pop))
rownames(richness.matrix) <- rownames(private.matrix) <- as.character(levels(test.data$loc.fac))

# Loop over each locus, doing richness and private calculations for each
for (locus.index in 1:length(levels(test.data$loc.fac))){
  # Extract the name of the current locus
  locus.id <- levels(test.data$loc.fac)[locus.index]
  # Subset data for that one locus
  locus.data <- as.data.frame(test.data$tab[, test.data$loc.fac == locus.id])
  # TODO: This needs only happen once?
  locus.data$pop <- test.data$pop
  
  # Do the allele counts for the locus
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
  
  # Transpose for appropriate allele x pop format needed by functions
  N.matrix <- t(count.matrix)
  richness.matrix[locus.index, ] <- calcRichness.all(N = N.matrix, g = 14)
  private.matrix[locus.index, ] <- calcPrivate.all(N = N.matrix, g = 14)
}

################################################################################
richness.test <- c(1, 1, 1, 1, 
                   1, 1, 1, 1, 
                   2, 1, 1, 1, 
                   1, 1, 1, 1, 
                   1, 1, 2, 1, 
                   1, 1, 1, 1, 
                   1, 1, 1, 1, 
                   2, 1, 1, 1, 
                   1, 1, 1, 1, 
                   1, 1, 1, 1)
private.test <- rep(x = 0, times = 40)

################################################################################
pass <- all(as.vector(richness.matrix) == richness.test) &
  all(as.vector(private.matrix) == private.test)

cat("Test 6 pass = ", pass, "\n", sep = "")

rm(list = ls())