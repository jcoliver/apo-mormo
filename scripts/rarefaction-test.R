# Rarefaction testing
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-05-17

################################################################################
rm(list = ls())

source(file = "functions/rarefaction.R")

################################################################################
## TEST 1
## two alleles, two populations
N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)

# Prob allele *not* found
calcQ.v(N.col = N.matrix[, 1], g = 2)
calcQ.v(N.col = N.matrix[, 1], g = 4)
calcQ.v(N.col = N.matrix[, 2], g = 2)
calcQ.v(N.col = N.matrix[, 2], g = 4)

# Prob allele found
calcP.v(N.col = N.matrix[, 1], g = 2)
calcP.v(N.col = N.matrix[, 1], g = 4)
calcP.v(N.col = N.matrix[, 2], g = 2)
calcP.v(N.col = N.matrix[, 2], g = 4)

# Prob for single allele **doesn't work if N.col = 0 (for allele not found in population)
# CANNOT work, because Nj (total # of genes surveyed from population is not included)
calcP.v(N.col = N.matrix[1, 2], g = 2)
calcP.v(N.col = N.matrix[2, 2], g = 2)

# Allele richness for individual populations
calcRichness(N.col = N.matrix[, 1], 4)
calcRichness(N.col = N.matrix[, 2], 4)

# Allele richness for each population
calcRichness.all(N = N.matrix, g = 2)
calcRichness.all(N = N.matrix, g = 4)

# Private alleles for each population
source(file = "functions/rarefaction.R")
calcPrivate(N = N.matrix, g = 2, j = 1)
calcPrivate(N = N.matrix, g = 4, j = 1)
calcPrivate(N = N.matrix, g = 2, j = 2)
calcPrivate(N = N.matrix, g = 4, j = 2)

# Private alleles for all populations
calcPrivate.all(N = N.matrix, g = 2)
calcPrivate.all(N = N.matrix, g = 4)

################################################################################
## TEST 2
## three alleles, four populations
source(file = "functions/rarefaction.R")

N.matrix <- matrix(data = c(3, 6, 4, 0, 
                            3, 0, 0, 0, 
                            0, 0, 0, 4), nrow = 3, ncol = 4, byrow = TRUE)
calcPrivate.all(N = N.matrix, g = 2)
calcPrivate.all(N = N.matrix, g = 4)

# Prob allele *not* found
calcQ.v(N.col = N.matrix[, 1], g = 2)
calcQ.v(N.col = N.matrix[, 1], g = 4)
calcQ.v(N.col = N.matrix[, 2], g = 4)
calcQ.v(N.col = N.matrix[, 3], g = 4)
calcQ.v(N.col = N.matrix[, 4], g = 4)

# Prob allele found
calcP.v(N.col = N.matrix[, 1], g = 2)
calcP.v(N.col = N.matrix[, 1], g = 4)
calcP.v(N.col = N.matrix[, 2], g = 2)
calcP.v(N.col = N.matrix[, 2], g = 4)

################################################################################
## TEST 3
## Work on genotype data
source(file = "functions/rarefaction.R")

# Can create N.matrix from using table on a data frame
# Example:
locus.data <- data.frame(pop = c(rep(x = 1, times = 6), rep(x = 2, times = 6), rep(x = 3, times = 4), rep(x = 4, times = 4)),
                         allele = c(rep(x = "10", times = 3), rep(x = "11", times = 3), rep(x = "10", times = 6), rep(x = "10", times = 4), rep(x = "12", times = 4)))
N.matrix <- t(unclass(table(locus.data)))
calcPrivate.all(N = N.matrix, g = 2)
calcPrivate.all(N = N.matrix, g = 4)
calcRichness.all(N = N.matrix, g = 2)
calcRichness.all(N = N.matrix, g = 4)

################################################################################
## TEST 3
## Work on real data, extracted from genind object
source(file = "functions/rarefaction.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

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

library("dplyr")

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
private <- calcPrivate.all(N = N.matrix, g = 14)
richness <- calcRichness.all(N = N.matrix, g = 14)
# Eventually give the rowname based on locus.id
# rownames(private) <- rownames(richness) <- locus.id

# Now want this for every locus, add row to matrix for each locus

################################################################################
## TEST 4
## Work on larger set of real data, extracted from genind object
source(file = "functions/rarefaction.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

# Grab some data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[1:35, 1:20]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:20])
test.data$pop <- factor(apo.str.genind@pop[1:35])

# Get allele x pop matrix for first locus
locus.data <- as.data.frame(test.data$tab[, 1:2])
locus.data$pop <- test.data$pop
locus.fac <- factor(test.data$loc.fac[1:2])
locus.id <- unique(test.data$loc.fac[1:2])

library("dplyr")

# Get counts for each column (allele) of the data
allele.counts <- locus.data %>% 
  group_by(pop) %>% 
  summarize_all(funs(sum))

# Pull out counts (first column is pop id)
count.matrix <- as.matrix(allele.counts[, c(2:ncol(allele.counts))])

# Set row names from values of pop (odd syntax to extract one column 
# from the allele.counts tibble)
rownames(count.matrix) <- as.character(allele.counts[[1]])
colnames(count.matrix) <- gsub(pattern = paste0(as.character(locus.id), "."),
                               replacement = "", 
                               x = colnames(count.matrix))

N.matrix <- t(count.matrix)
private <- calcPrivate.all(N = N.matrix, g = 14)
richness <- calcRichness.all(N = N.matrix, g = 14)
# Eventually give the rowname based on locus.id
# rownames(private) <- rownames(richness) <- locus.id

################################################################################
## TEST 5
## Subset of real data, looping over all loci

source(file = "functions/rarefaction.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

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
## TEST 6
## Subset of loci for all individuals, looping over all loci

source(file = "functions/rarefaction.R")
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

# Subset data for testing (using a list not an S4 object)
test.data = list()
test.data$tab <- apo.str.genind@tab[, 1:20]
test.data$loc.fac <- factor(apo.str.genind@loc.fac[1:20])
test.data$pop <- factor(apo.str.genind@pop)

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
## TEST 7
## Investigating troublesome populations 9, 11, 12, 13

################################################################################
## TEST 8
## Subset of loci for all individuals, looping over all loci using function

rarefied.data <- rarefiedMatrices(data = test.data, g = 6)

#' Now we need a means of getting that kind of data frame extracted from genind
#' object. Ultimately want: alleles x pop
#'        pop
#' allele 1 2 3 4
#'     10 3 6 4 0
#'     11 3 0 0 0
#'     12 0 0 0 4
#' The two things to look at are:
#' genind@tab: a table of individual genotype counts
#' apo.str.genind@tab[1:10, 1:6]
#            L0001.10 L0001.12 L0002.13 L0002.11 L0003.13 L0003.11
# HullMtn_48        2        0        2        0        2        0
# HullMtn_50        2        0        2        0        2        0
# HullMtn_51        2        0        2        0        2        0
# HullMtn_52        2        0        2        0        2        0
# HullMtn_53        2        0        2        0        2        0
# HullMtn_54        2        0        2        0        2        0
# HullMtn_55        2        0        2        0        2        0
# Ladoga_82         2        0        2        0        2        0
# Ladoga_83         2        0        2        0        2        0
# Ladoga_85         2        0        2        0        2        0
#' The locus values can be seen via genind@loc.fac:
#' L0001 L0001 L0002 L0002 L0003 L0003
#' These are indexed in the same order as columns appear in genind@tab
#' genind@pop: a vector of type factor indicating population membership; it is 
#' indexed identically to the rows of genind@tab
#' apo.str.genind@pop[1:10]
# 1  2  3  4  5  6  7  8  9 10 
# 1  1  1  1  1  1  1  2  2  2 
# Levels: 1 11 12 13 14 2 3 4 5 6 7 8 9

