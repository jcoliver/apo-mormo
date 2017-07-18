# Permutation tests of rarefied allele richness and private allele counts
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-18

rm(list = ls())

################################################################################
# Load dependencies
library("adegenet")
source(file = "functions/rarefaction.R")

# Load data
genind.file = "output/genind-object.RData"
load(file = genind.file) # apo.str.genind

################################################################################
# Extract the data to use
data.list = list()
data.list$tab <- apo.str.genind@tab
data.list$loc.fac <- factor(apo.str.genind@loc.fac)
data.list$pop <- factor(apo.str.genind@pop)

################################################################################
# Run rarefaction analyses, using smallest population as minimal cutoff
min.size <- min(table(data.list$pop))
rarefied.values <- rarefiedMatrices(data = data.list, 
                                    g = min.size, 
                                    display.progress = TRUE)
richness.matrix <- rarefied.values$richness
private.matrix <- rarefied.values$private

################################################################################
# Calculate sums
richness.sums <- apply(X = richness.matrix, 
                       MARGIN = 2, 
                       FUN = sum, 
                       na.rm = TRUE)

private.sums <- apply(X = private.matrix,
                      MARGIN = 2, 
                      FUN = sum, 
                      na.rm = TRUE)

################################################################################
# Loop to create distribution of sums
num.perm <- 100
print.freq <- 10

richness.perms <- matrix(data = 0, 
                         nrow = num.perm, 
                         ncol = length(levels(data.list$pop)))
private.perms <- matrix(data = 0, 
                        nrow = num.perm, 
                        ncol = length(levels(data.list$pop)))

colnames(richness.perms) <- colnames(private.perms) <- levels(data.list$pop)
rownames(richness.perms) <- rownames(private.perms) <- 1:num.perm

for (perm in 1:num.perm) {
  if (perm %% print.freq == 0 || perm == 1 || perm == num.perm) {
    cat("Running permutation ", perm, "\n", sep = "")
  }
  shuffle.data <- list()
  shuffle.data <- data.list
  shuffle.data$tab <- shuffle.data$tab[sample(nrow(shuffle.data$tab)), ]
  rownames(shuffle.data$tab) <- rownames(data.list$tab)
  
  shuffle.min.size <- min(table(shuffle.data$pop))
  shuffle.rarefied.values <- rarefiedMatrices(data = shuffle.data, 
                                              g = min.size, 
                                              display.progress = FALSE)
  shuffle.richness.matrix <- shuffle.rarefied.values$richness
  shuffle.private.matrix <- shuffle.rarefied.values$private
  
  richness.perms[perm, ] <- apply(X = shuffle.richness.matrix,
                                 MARGIN = 2,
                                 FUN = sum,
                                 na.rm = TRUE)
  private.perms[perm, ] <- apply(X = shuffle.private.matrix,
                                MARGIN = 2,
                                FUN = sum,
                                na.rm = TRUE)
}

################################################################################
# From permutation distributions assess significance of observed values
richness.sig <- data.frame(direction = rep(x = NA, times = length(levels(data.list$pop))),
                           p.value = 0)
richness.sig$direction <- factor(richness.sig$direction, levels = c("High", "Low"))
rownames(richness.sig) <- levels(data.list$pop)
private.sig <- richness.sig

for (one.pop in 1:ncol(richness.perms)) {
  richness.perm.result <- permSignificance(obs.data = richness.sums[one.pop], richness.perms[, one.pop])
  richness.sig$direction[one.pop] <- richness.perm.result$direction
  richness.sig$p.value[one.pop] <- richness.perm.result$p.value
  private.perm.result <- permSignificance(obs.data = private.sums[one.pop], private.perms[, one.pop])
  private.sig$direction[one.pop] <- private.perm.result$direction
  private.sig$p.value[one.pop] <- private.perm.result$p.value
}

richness.sig
private.sig
