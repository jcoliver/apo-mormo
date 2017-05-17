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
calcQ <- function(Nj, Nij, g) {
  if (Nj <= (Nij + g - 1)) {
    return(0)
  }
  prod <- 1
  iter.g <- 0
  while (iter.g < g && prod > 0) {
    numerator <- Nj - Nij - iter.g
    if (numerator == 0) {
      prod <- 0
    } else {
      denominator <- Nj - iter.g
      prod <- prod * (numerator / denominator)
      iter.g <- iter.g + 1
    }
  }
  return(prod)
}

################################################################################
# Probability of observing allele i in population j given sample of g genes
# Nj <- sum(N[ ,j])
# Nij <- N[i, j]
calcP <- function(Nj, Nij, g) {
  prob <- 1 - calcQ(Nj = Nj, Nij = Nij, g = g)
  return(prob)
}

################################################################################
# Probability of *not* observing allele in population given sample of g 
# genes
# N.col is vector of allele counts for population j
calcQ.v <- function(N.col, g) {
  if(length(N.col) < 2 || sum(N.col) == 0) {
    return(NA)
  }
  Nj <- sum(N.col)
  product <- rep(x = 1, times = length(N.col))
  iter.g <- 0
  while (iter.g < g && sum(product) > 0) {
    numerator <- Nj - N.col - iter.g
    denominator <- Nj - iter.g
    product <- product * (numerator / denominator)
    iter.g <- iter.g + 1
  }
  return(product)
}

################################################################################
# Probability of observing allele in population given sample of g genes
# N.col is vector of allele counts for population j
calcP.v <- function(N.col, g) {
  prob <- 1 - calcQ.v(N.col = N.col, g = g)
  return(prob)
}

printArgs <- function(Nj, Nij, g, standin.for.nij) {
  Nij <- standin.for.nij
  cat("Nj = ", Nj, "\n")
  cat("Nij = ", Nij, "\n")
  cat("g = ", g, "\n")
  return(calcP(Nj = Nj, Nij = Nij, g = g))
}

################################################################################
# Make N be a vector, not a matrix
# calcRichness <- function(N, j, g) {
calcRichness <- function(N, j, g) {
  if(require(package = "plyr")) {
    stop("ERROR: Richness calculations requires plyr package.")
  }
  richness <- 0
  m <- nrow(N)
  for (i in 1:m) {
    richness <- richness + calcP(Nj = sum(N[ ,j]), Nij = N[i, j], g = g)
  }
  # return(richness)

  aaply(.data = N, .margins = 2, .fun = calcP, Nj = 1, Nij = 2, g = 2)
  
  ## TODO: here, but stumbling on apply / aaply
  ## Try first with for loop, then attempt a *ply approach
  # apply(X = N, MARGIN = 1, FUN = printArgs, Nj = 1, Nij = 2, g = 2)
  # richness <- aaply(.data = N, .margins = 2, .fun = printArgs, Nj = sum(N[, j]), g = 2, standin.for.nij = N[i, j])
  # return(sum(richness))
  # p.vector <- apply(X = N, MARGIN = 1, FUN = calcP, Nj = sum(x), Nij = x[j], g = g)
  # return(sum(p.vector))
}

