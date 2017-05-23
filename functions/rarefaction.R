# Private allele rarefaction develoment
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-05-16

################################################################################
# SUMMARY
# Kalinowski's algorithms for correcting for uneven sampling to estimate allele 
# richness and private allele frequency

################################################################################
#' Probability of not observing alleles in single population given sample of g 
#'   genes
#' @param N.col vector of allele counts for a single population
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates the probability of not observing each allele in a 
#'   given population for a rarefied sample of size \code{g}
#' @return vector of probabilities for single population
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
#' Rarefied allele probability
#' @param N.col vector of allele counts for a single population
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates the probability of observing each allele in a given 
#'   population for a rarefied sample of size \code{g}
#' @return vector of probabilities for single population
calcP.v <- function(N.col, g) {
  prob <- 1 - calcQ.v(N.col = N.col, g = g)
  return(prob)
}

################################################################################
#' Rarefied allele richness for a single population
#' @param N.col vector of allele counts for a single population
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates allele richness for a rarefied sample of size 
#'   \code{g}
#' @return Allele richness in a single population rarefied to size \code{g}
calcRichness <- function(N.col, g) {
  sum.richness <- sum(calcP.v(N.col = N.col, g = g))
  return(sum.richness)
}

################################################################################
#' Rarefied allele richness for each population
#' @param N matrix of allele counts for all populations, with alleles as rows 
#'   and populations as columns
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates allele richness for a rarefied sample of size 
#'   \code{g}
#' @return A vector of allele richness, where each element is the richness for 
#'   a single population indexed by columns in \code{N}
calcRichness.all <- function(N, g) {
  all.rich <- apply(X = N, MARGIN = 2, FUN = function(x) {calcRichness(N.col = x, g = g)})
  return(all.rich)
}

################################################################################
#' Rarefied private allele count for single population
#' @param N matrix of allele counts for all populations, with alleles as rows 
#'   and populations as columns
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @param j index of the population for which to perform calculations
#' 
calcPrivate <- function(N, g, j) {
  # sum over all alleles
  # Pijg * [(]prod from j'=1 to J, j' != j (Qij'g)]
  m <- nrow(N)
  private.sum <- 0
  Pijg <- calcP.v(N.col = N[, j], g = g) # a vector of allele probabilities
  
  # Now we need a vector of the Q products
  Q.matrix <- apply(X = N, MARGIN = 2, FUN = function(x) {calcQ.v(N.col = x, g = g)})
  Q.matrix[, j] <- 1 # Dummy coding, so product of row is unaffected for column j
  Q.products <- apply(X = Q.matrix, MARGIN = 1, FUN = prod)

  # Have vector Pijg and vector of product Qij'g
  # Multiply corresponding elements in each, then sum result
  pi.hat <- sum(Pijg * Q.products)
  return(pi.hat)
}
