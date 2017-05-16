# Private allele rarefaction develoment
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-05-16

################################################################################
# SUMMARY
# Use Kalinowski's algorithm for correcting for uneven sampling to estimate
# private allele frequency

# calcQ <- function(N, i, j, g) {
#   Nj <- sum(N[ ,j])
#   Nij <- N[i, j]
#   prod <- 1
#   iter.g <- 0
#   while (iter.g < g && prod > 0) {
#     numerator <- Nj - Nij - iter.g
#     if (numerator == 0) {
#       prod <- 0
#     } else {
#       denominator <- Nj - iter.g
#       prod <- prod * (numerator / denominator)
#       iter.g <- iter.g + 1
#     }
#   }
#   return(prod)
# }

# calcP <- function(N, i, j, g) {
#   prob <- 1 - calcQ(N = N, i = i, j = j, g = g)
#   return(prob)
# }

# # P112 = 1 - Q112
# calcP(N = N.matrix, i = 1, j = 1, g = 2)
# 
# # Q122
# calcQ(N = N.matrix, i = 1, j = 2, g = 2)
# 
# # P212 = 1 - Q212
# calcP(N = N.matrix, i = 2, j = 1, g = 2)
# 
# # Q222
# calcQ(N = N.matrix, i = 2, j = 2, g = 2)

calcQ <- function(Nj, Nij, g) {
  # Nj <- sum(N[ ,j])
  # Nij <- N[i, j]
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

calcP <- function(Nj, Nij, g) {
  prob <- 1 - calcQ(Nj = Nj, Nij = Nij, g = g)
  return(prob)
}

calcRichness <- function(N, j, g) {
  p.vector <- apply(X = N, MARGIN = 1, FUN = calcP, i = 1, j = j, g = g)
  return(sum(p.vector))
}

N.matrix <- matrix(data = c(3, 6, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)

# P112 = 1 - Q112
calcP(Nj = sum(N.matrix[, 1]), Nij = N.matrix[1, 1], g = 2)

# Q122
calcQ(Nj = sum(N.matrix[, 2]), Nij = N.matrix[1, 2], g = 2)

# P212 = 1 - Q212
calcP(Nj = sum(N.matrix[, 1]), Nij = N.matrix[2, 1], g = 2)

# Q222
calcQ(Nj = sum(N.matrix[, 2]), Nij = N.matrix[2, 2], g = 2)
