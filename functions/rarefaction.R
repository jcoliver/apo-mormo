# Private allele rarefaction develoment
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-05-16

################################################################################
# SUMMARY
# Kalinowski's algorithms for correcting for uneven sampling to estimate allele 
# richness and private allele frequency

################################################################################
#' Probability of allele absence in rarefied sample
#' @param N.col vector of allele counts for a single population
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates the probability of not observing each allele in a 
#'   given population for a rarefied sample of size \code{g}
#' @return vector of probabilities for single population
#' @examples 
#' N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' # Allele probabilities for population 1 given sample of 2 genes
#' calcQ.v(N.col = N.matrix[, 1], g = 2)
#' # Allele probabilities for population 2 given sample of 4 genes
#' calcQ.v(N.col = N.matrix[, 2], g = 4)
calcQ.v <- function(N.col, g) {
  if(length(N.col) < 2 || sum(N.col) == 0 || sum(N.col) < g) {
    return(rep(x = NA, times = length(N.col)))
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
#' Probability of allele presence in rarefied sample
#' @param N.col vector of allele counts for a single population
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates the probability of observing each allele in a given 
#'   population for a rarefied sample of size \code{g}
#' @return vector of probabilities for single population
#' @examples 
#' N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' # Allele probabilities for population 1 given sample of 2 genes
#' calcP.v(N.col = N.matrix[, 1], g = 2)
#' # Allele probabilities for population 2 given sample of 4 genes
#' calcP.v(N.col = N.matrix[, 2], g = 4)
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
#' @examples 
#' N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' # Allele richness for population 1 given sample of 2 genes
#' calcRichness.all(N.col = N.matrix[, 1], g = 2)
#' # Allele richness for population 2 given sample of 4 genes
#' calcRichness.all(N.col = N.matrix[, 2], g = 4)
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
#' @examples 
#' N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' # Allele richness for all populations given sample of 4 genes
#' calcRichness.all(N = N.matrix, g = 4)
calcRichness.all <- function(N, g) {
  all.rich <- apply(X = N, MARGIN = 2, FUN = function(x) {calcRichness(N.col = x, g = g)})
  # Preserve population names if they exist in N
  names(all.rich) <- colnames(N)
  return(all.rich)
}

################################################################################
#' Rarefied private allele count for single population
#' @param N matrix of allele counts for all populations, with alleles as rows 
#'   and populations as columns
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @param j index of the population for which to perform calculations
#' @description Calculates the expected number of private alleles in population
#' \code{j} given a sample of \code{g} genes
#' @return A vector of length 1, with the expected number of private alleles for
#' the given population
#' @examples 
#' N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' # Private alleles in population 1, for a sample size of 2 genes
#' calcPrivate(N = N.matrix, g = 2, j = 1)
#' # Private alleles in population 2, for a sample size of 4 genes
#' calcPrivate(N = N.matrix, g = 4, j = 2)
calcPrivate <- function(N, g, j) {
  # TODO: update to mirror calcPrivate.all treatment of populations with 
  # sampling below g; will affect *all* populations if any one population
  # has fewer than g genes sampled
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

################################################################################
#' Rarefied private allele count for all populations
#' @param N matrix of allele counts for all populations, with alleles as rows 
#'   and populations as columns
#' @param g gene sample size; generally the number of haploids or twice the 
#'   number of diploid individuals in the smallest population sample size
#' @description Calculates the expected number of private alleles for a sample 
#' of \code{g} genes
#' @return A vector of \code{nrow(N)}, with the expected number of private 
#' alleles for each population
#' @examples 
#' N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' # Private alleles for sample of 2 genes
#' calcPrivate.all(N = N.matrix, g = 2)
calcPrivate.all <- function(N, g) {
  # Start by reducing matrix to only those populations with at least g sampled 
  # genes
  cols.keep <- which(colSums(x = N, na.rm = TRUE) >= g)
  N.analyze <- N[, cols.keep]
  private.alleles.return <- rep(NA, times = ncol(N))
  
  if (is.matrix(N.analyze) && dim(N.analyze)[1] > 1) { # Need at least two populations for this calculation?
    # Both P and Q matrix have alleles in rows, and populations in columns
    #' P114 P124 P134 P144
    #' P214 P224 P234 P244
    #' P314 P324 P334 P344
    #' P414 P424 P434 P444
    P.matrix <- apply(X = N.analyze, MARGIN = 2, FUN = function(x) {calcP.v(N.col = x, g = g)})
    Q.matrix <- apply(X = N.analyze, MARGIN = 2, FUN = function(x) {calcQ.v(N.col = x, g = g)})
    
    # Private allele calculation:
    #' For population 1, it is the sum of (three alleles, four populations):
    #' P114 x (Q124 x Q134 x Q144)
    #' P214 x (Q224 x Q234 x Q244)
    #' P314 x (Q324 x Q334 x Q344)
    #' 
    #' For population 2, it is the sum of (three alleles, four populations):
    #' P124 x (Q114 x Q134 x Q144)
    #' P224 x (Q214 x Q234 x Q244)
    #' P324 x (Q314 x Q334 x Q344)
    
    # private allele count vector; each element corresponds to a population
    private.alleles <- numeric(ncol(N.analyze))
    for (j in 1:ncol(N.analyze)) {
      corrected.Q.matrix <- Q.matrix
      corrected.Q.matrix[, j] <- 1 # Dummy coding, so product of row is unaffected for column j
      Q.products <- apply(X = corrected.Q.matrix, MARGIN = 1, FUN = prod)
      private.alleles[j] <- sum(P.matrix[, j] * Q.products)
    }
    
    # Now we need a vector of length equivalent to the original number of columns 
    # in N, including those with too few sampled genes
    private.alleles.return[cols.keep] <- private.alleles
    
  }  
  # Preserve population names if they exist in N
  names(private.alleles.return) <- colnames(x = N)
  return(private.alleles.return)
}


################################################################################
#' Private allele and allele richness for all loci and all populations
#' 
rarefiedMatrices <- function(data, g = 2, display.progress = FALSE) {
  if (!require("dplyr")) {
    stop("rarefiedMatrices requires the dplyr package.")
  }
  
  # Establish matrices that will hold the final results
  richness.matrix <- matrix(data = 0, 
                            nrow = length(levels(data$loc.fac)),
                            ncol = length(levels(data$pop)))
  
  private.matrix <- matrix(data = 0, 
                           nrow = length(levels(data$loc.fac)),
                           ncol = length(levels(data$pop)))
  
  colnames(richness.matrix) <- colnames(private.matrix) <- as.character(levels(data$pop))
  rownames(richness.matrix) <- rownames(private.matrix) <- as.character(levels(data$loc.fac))
  
  if (display.progress) {
    progress.bar <- txtProgressBar(min = 1, max = length(levels(data$loc.fac)), style = 3)
  }
  # Loop over each locus, doing richness and private calculations for each
  for (locus.index in 1:length(levels(data$loc.fac))){
    if (display.progress) {
      setTxtProgressBar(pb = progress.bar, value = locus.index)
    }
    
    # Extract the name of the current locus
    locus.id <- levels(data$loc.fac)[locus.index]
    # Subset data for that one locus
    locus.data <- as.data.frame(data$tab[, data$loc.fac == locus.id])
    # TODO: This needs only happen once?
    locus.data$pop <- data$pop
    
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
    richness.matrix[locus.index, ] <- calcRichness.all(N = N.matrix, g = g)
    private.matrix[locus.index, ] <- calcPrivate.all(N = N.matrix, g = g)
  }
  close(progress.bar)
  return(list(richness = richness.matrix, private = private.matrix))
}