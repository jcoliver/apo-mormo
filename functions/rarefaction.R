# Private allele rarefaction develoment
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-05-16

################################################################################
# SUMMARY
# Kalinowski's algorithms for correcting for uneven sampling to estimate allele 
# richness and private allele frequency

# Probability of *not* observing allele i in population j given sample of g 
# genes
# Nj <- sum(N[ ,j])
# Nij <- N[i, j]
# DEPRECATED
calcQ <- function(Nj, Nij, g) {
  if (Nj <= (Nij + g - 1)) {
    return(0)
  }
  prod <- 1
  u <- 0
  while (u < g && prod > 0) {
    numerator <- Nj - Nij - u
    if (numerator == 0) {
      prod <- 0
    } else {
      denominator <- Nj - u
      prod <- prod * (numerator / denominator)
      u <- u + 1
    }
  }
  return(prod)
}

################################################################################
# Probability of observing allele i in population j given sample of g genes
# Nj <- sum(N[ ,j])
# Nij <- N[i, j]
# DEPRECATED
calcP <- function(Nj, Nij, g) {
  prob <- 1 - calcQ(Nj = Nj, Nij = Nij, g = g)
  return(prob)
}

################################################################################
# Calculate allele richness for population j given a sample of g genes
# N is matrix of allele counts for all populations, with alleles as rows and 
# populations as columns
# DEPRECATED
calcRichness <- function(N, j, g) {
  richness <- 0
  m <- nrow(N)
  for (i in 1:m) {
    richness <- richness + calcP(Nj = sum(N[ ,j]), Nij = N[i, j], g = g)
  }
  return(richness)
}

################################################################################
# Probability of *not* observing allele in single population given sample of g 
# genes
# N.col is vector of allele counts for population j
calcQ.v <- function(N.col, g) {
  if(length(N.col) < 2 || sum(N.col) == 0) {
    return(NA)
  }
  Nj <- sum(N.col)
  product <- rep(x = 1, times = length(N.col))
  u <- 0

  # There is a quicker way to do this, given:
  #' For each allele, the numerator is an integer range: 
  #' (Nj - Nij) to (Nj - Nij - (g - 1))
  #' And the corresponding denominator range (identical across alleles):
  #' (Nj) to (Nj - (g - 1))
  #' It is easy to create the vector of denominators:
  #'   denominators <- seq(from = Nj, to = Nj - (g - 1))
  #' But numerators will need to be a matrix...
  
  while (u < g && sum(product) > 0) {
    numerator <- Nj - N.col - u
    denominator <- Nj - u
    product <- product * (numerator / denominator)
    u <- u + 1
  }
  return(product)
}

################################################################################
# Probability of observing allele in single population given sample of g genes
# N.col is vector of allele counts for population j
calcP.v <- function(N.col, g) {
  prob <- 1 - calcQ.v(N.col = N.col, g = g)
  return(prob)
}

################################################################################
# Calculate allele richness for a single population given sample of g genes
# N.col is vector of allele counts for population j
calcRichness.v <- function(N.col, g) {
  sum.richness <- sum(calcP.v(N.col = N.col, g = g))
  return(sum.richness)
}

################################################################################
# Calculate allele richness in each population
# N is matrix of allele counts for all populations, with alleles as rows and 
# populations as columns
#' Calculate allele richness for each population
#' @param N matrix of allele counts for all populations, with alleles as rows 
#' and populations as columns
#' @param g gene sample size; generally the number of haploids or twice the 
#' number of diploid individuals in the smallest population sample size
#' @description Calculates allele richness for a rarefied sample of size 
#' \code{g}
#' @return A vector of allele richness, where each element is the richness for 
#' a single population indexed by columns in \code{N}
calcRichness.all <- function(N, g) {
  all.rich <- apply(X = N, MARGIN = 2, FUN = function(x) {calcRichness.v(N.col = x, g = g)})
  return(all.rich)
}
