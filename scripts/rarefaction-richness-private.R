# Allele richness and private allele counts using rarefaction
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-07-18

rm(list = ls())

library("adegenet")

source(file = "functions/rarefaction-functions.R")
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
min.size <- min(table(data.list$pop)) * 2
rarefied.values <- rarefiedMatrices(data = data.list, 
                                    g = min.size, 
                                    display.progress = TRUE)
richness.matrix <- rarefied.values$richness
private.matrix <- rarefied.values$private

################################################################################
# Extract population means
richness.means <- apply(X = richness.matrix, 
                        MARGIN = 2, 
                        FUN = mean, 
                        na.rm = TRUE)

private.means <- apply(X = private.matrix,
                      MARGIN = 2, 
                      FUN = mean, 
                      na.rm = TRUE)

################################################################################
# Extract population sd
richness.sd <- apply(X = richness.matrix, 
                     MARGIN = 2, 
                     FUN = sd, 
                     na.rm = TRUE)

private.sd <- apply(X = private.matrix,
                    MARGIN = 2, 
                    FUN = sd, 
                    na.rm = TRUE)

################################################################################
# Export to file
write.table(x = data.frame(population = names(richness.means), 
                           richness.means,
                           richness.sd, 
                           private.means, 
                           private.sd),
            file = "output/rarefied-alleles-observed.csv",
            sep = ",",
            row.names = FALSE,
            quote = FALSE)
