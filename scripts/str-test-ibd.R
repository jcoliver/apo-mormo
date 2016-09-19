# Testing import of STRUCTURE file
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-07

library("adegenet")
library("hierfstat")

################################################################################
# Current structure file, data/Apodemia_0.9-noDockweiler-GNPSK.str
apo.str.genind <- read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str",
                                 n.ind = 102,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")
#' /// GENIND OBJECT /////////
#'   
#'   // 102 individuals; 4,057 loci; 7,808 alleles; size: 4.8 Mb
#' 
#' // Basic content
#' @tab:  102 x 7808 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 1-2)
#' @loc.fac: locus factor for the 7808 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str", 
#'                       n.ind = 102, n.loc = 4057, onerowperind = TRUE, col.lab = 1, 
#'                       col.pop = 2, col.others = 0, row.marknames = 0, NA.char = "-9", 
#'                       sep = "\t")
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 4-12)
#' @other: a list containing: X 

apo.str.summary <- summary(apo.str.genind)
mean(apo.str.summary$Hexp)
# 0.097
mean(apo.str.summary$Hobs)
# 0.051

start <- Sys.time()
curr.str.p.fst <- genet.dist(apo.str.genind, method = "WC84") # ~10.5 minutes
end <- Sys.time()
max(p.fst)

min(p.fst)

mean(p.fst)


################################################################################
# Much smaller, earlier structure file data/filtered_Apo.str
apo.str.old.genind <- read.structure(file = "data/filtered_Apo.str",
                                 n.ind = 109,
                                 n.loc = 432,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")

#' /// GENIND OBJECT /////////
#'   
#'   // 109 individuals; 432 loci; 864 alleles; size: 585.7 Kb
#' 
#' // Basic content
#' @tab:  109 x 864 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 2-2)
#' @loc.fac: locus factor for the 864 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: read.structure(file = "data/filtered_Apo.str", n.ind = 109, n.loc = 432, 
#'                       onerowperind = TRUE, col.lab = 1, col.pop = 2, col.others = 0, 
#'                       row.marknames = 0, NA.char = "-9", sep = "\t")
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 2-12)
#' @other: a list containing: X 

apo.str.old.summary <- summary(apo.str.old.genind)
mean(apo.str.old.summary$Hexp)
# 0.101
mean(apo.str.old.summary$Hobs)
# 0.052

start <- Sys.time()
old.str.p.fst <- genet.dist(apo.str.old.genind, method = "WC84") # ~1.3 minutes
end <- Sys.time()
max(old.str.p.fst)
# 0.81
min(old.str.p.fst)
# 0.02
mean(old.str.p.fst)
# 0.45

old.nei.fst <- genet.dist(dat = apo.str.old.genind, method = "Nei87")

################################################################################
# Current structure file plus Dockweiler & Grasslands NPSK populations, 
# data/Apodemia_0.9-noDockweiler-GNPSK.str
apo.str.dg.genind <- read.structure(file = "data/Apodemia_0.9miss.str",
                                 n.ind = 109,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")
#' /// GENIND OBJECT /////////
#'   
#'   // 102 individuals; 4,057 loci; 7,808 alleles; size: 4.8 Mb
#' 
#' // Basic content
#' @tab:  102 x 7808 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 1-2)
#' @loc.fac: locus factor for the 7808 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str", 
#'                       n.ind = 102, n.loc = 4057, onerowperind = TRUE, col.lab = 1, 
#'                       col.pop = 2, col.others = 0, row.marknames = 0, NA.char = "-9", 
#'                       sep = "\t")
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 4-12)
#' @other: a list containing: X 

apo.str.dg.summary <- summary(apo.str.dg.genind)
mean(apo.str.dg.summary$Hexp)
# 0.100
mean(apo.str.dg.summary$Hobs)
# 0.049

start <- Sys.time()
curr.dg.str.p.fst <- genet.dist(apo.str.dg.genind, method = "WC84") # ~10.5 minutes
end <- Sys.time()
end - start
max(curr.dg.str.p.fst)
# 0.796
min(curr.dg.str.p.fst)
# 0.049
mean(curr.dg.str.p.fst)
# 0.461

################################################################################
apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     header = FALSE, 
                     stringsAsFactors = FALSE)

# 20 loci
source(file = "functions/apodemia-functions.R")
compare.20 <- compare.df.str(data = subset, inds = c(1:nrow(apo.df)), n.loci = 20)
  



#inds <- c(1:3, 80:82, 91:93)
inds <- c(1:20) # HullMtn, Lagoda, & some NEVallejo
n.loci <- 10

small.df <- apo.df[inds, 1:(2 + n.loci * 2)]

write.table(x = small.df, 
            file = "data/Apodemia_0.9-small-test.str",
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE,
            col.names = FALSE)

library("adegenet")
apo.str.genid <- read.structure(file = "data/Apodemia_0.9-small-test.str",
                              n.ind = length(inds),
                              n.loc = n.loci,
                              onerowperind = TRUE,
                              col.lab = 1,
                              col.pop = 2,
                              col.others = 0,
                              row.marknames = 0,
                              NA.char = "-9",
                              sep = "\t")

# Pairwise Fst
library("hierfstat")
p.fst.str <- genet.dist(apo.str.genid, method="WC84") # Too high?

# Combine columns for adegenet-readable df
#for.adegenet <- data.frame(inds = small.df[, 1], pops = small.df[, 2])
for.adegenet <- data.frame(small.df[, 0])
odd.cols <- seq(from = 3, to = (ncol(small.df) - 1), by = 2)
for (i in 1:length(odd.cols)) {
  odd.col <- odd.cols[i]
  even.col <- odd.col + 1
  genotype <- paste0(small.df[, odd.col], small.df[, even.col])
  locus = paste0("locus", i)
  genotype <- gsub(pattern = "-9-9", replacement = NA, x = genotype)
  for.adegenet[, locus] <- genotype
}
apo.df.genid <- df2genind(X = for.adegenet, 
                          ncode = 2, 
                          ind.names = small.df[, 1], 
                          pop = small.df[, 2], 
                          ploidy = 2)

p.fst.df <- genet.dist(apo.df.genid, method="WC84")

################################################################################
# Test on medium sized data set
apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     header = FALSE, 
                     stringsAsFactors = FALSE)


inds <- c(1:51) # half the individuals
n.loci <- floor(x = (ncol(apo.df) - 2) / 4) # ~ half the loci

medium.df <- apo.df[inds, 1:(2 + n.loci * 2)]

write.table(x = medium.df, 
            file = "data/Apodemia_0.9-medium-test.str",
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE,
            col.names = FALSE)

library("adegenet")
apo.str.genid <- read.structure(file = "data/Apodemia_0.9-medium-test.str",
                                n.ind = length(inds),
                                n.loc = n.loci,
                                onerowperind = TRUE,
                                col.lab = 1,
                                col.pop = 2,
                                col.others = 0,
                                row.marknames = 0,
                                NA.char = "-9",
                                sep = "\t")

# Pairwise Fst
library("hierfstat")
start <- Sys.time()
p.fst.str <- genet.dist(apo.str.genid, method="WC84") # ~1 minute
end <- Sys.time()
diff <- difftime(time1 = end, time2 = start, units = "auto")
diff
cat("Min: ", min(p.fst.str), ", Max: ", max(p.fst.str), "\n", sep = "")


# reading in from data.frame
for.adegenet <- data.frame(medium.df[, 0])
odd.cols <- seq(from = 3, to = (ncol(medium.df) - 1), by = 2)
for (i in 1:length(odd.cols)) {
  odd.col <- odd.cols[i]
  even.col <- odd.col + 1
  genotype <- paste0(medium.df[, odd.col], medium.df[, even.col])
  locus = paste0("locus", i)
  genotype <- gsub(pattern = "-9-9", replacement = NA, x = genotype)
  for.adegenet[, locus] <- genotype
}
apo.df.genid <- df2genind(X = for.adegenet, 
                          ncode = 2, 
                          ind.names = medium.df[, 1], 
                          pop = medium.df[, 2], 
                          ploidy = 2)

start <- Sys.time()
p.fst.df <- genet.dist(apo.df.genid, method="WC84") # ~ 1 minute
end <- Sys.time()
diff <- difftime(time1 = end, time2 = start, units = "auto")
diff
cat("Min: ", min(p.fst.df), ", Max: ", max(p.fst.df), "\n", sep = "")



################################################################################
# Test on the whole data set
apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     header = FALSE, 
                     stringsAsFactors = FALSE)


inds <- c(1:nrow(apo.df))
n.loci <- (ncol(apo.df) - 2) / 2

full.df <- apo.df[inds, 1:(2 + n.loci * 2)]

write.table(x = full.df, 
            file = "data/Apodemia_0.9-full-test.str",
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE,
            col.names = FALSE)

library("adegenet")
apo.str.genid <- read.structure(file = "data/Apodemia_0.9-full-test.str",
                                n.ind = length(inds),
                                n.loc = n.loci,
                                onerowperind = TRUE,
                                col.lab = 1,
                                col.pop = 2,
                                col.others = 0,
                                row.marknames = 0,
                                NA.char = "-9",
                                sep = "\t")

# Pairwise Fst
library("hierfstat")
start <- Sys.time()
p.fst.str <- genet.dist(apo.str.genid, method="WC84") # Takes ~10 minutes
end <- Sys.time()
diff <- difftime(time1 = end, time2 = start, units = "auto")
diff

################################################################################
# Compare 2-character vs. 1-character alleles (small)

apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     header = FALSE, 
                     stringsAsFactors = FALSE)

#inds <- c(1:15) # HullMtn & Lagoda
#loci <- c(13:22)
# inds <- c(1:23) # HullMtn, Lagoda, NEVallejo
# loci <- c(3:4000)
inds <- c(1:102) # all individuals
loci <- c(3:4060) # ~ half the loci

small.df <- apo.df[inds, loci]

# Prepare a 1-character matrix
small.df.1c <- small.df
# replace 10 with 0, 11 with 1...
replace.df <- data.frame(from = c(10, 11, 12, 13), to = c(0, 1, 2, 3))
for (i in 1:ncol(small.df.1c)) {
  for (j in 1:nrow(replace.df)) {
    matching.cells <- which(small.df.1c[, i] == replace.df[j, 1])
    small.df.1c[matching.cells, i] <- replace.df[j, 2]
  }
}

library("adegenet")
genid.2c <- df2genind(X = small.df,
                      ncode = 2,
                      ind.names = apo.df[inds, 1],
                      pop = apo.df[inds, 2],
                      ploidy = 2)

genid.1c <- df2genind(X = small.df.1c,
                      ncode = 1,
                      ind.names = apo.df[inds, 1],
                      pop = apo.df[inds, 2],
                      ploidy = 2)

library("hierfstat")
p.fst.2c <- genet.dist(genid.2c, method = "WC84")
p.fst.1c <- genet.dist(genid.1c, method = "WC84")
p.fst.2c
p.fst.1c

# To compare with full data set:
# source(file = "functions/apodemia-functions.R")
# p.str <- p.fst.str(filename = "data/Apodemia_0.9-noDockweiler-GNPSK.str", 
#                    n.ind = length(inds), 
#                    n.loci = length(loci))

################################################################################
# Create data frame were genotype (two alleles) is in a single column

apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     header = FALSE, 
                     stringsAsFactors = FALSE)

# inds <- c(1:15) # HullMtn & Lagoda
# loci <- c(13:22)
# inds <- c(1:23) # HullMtn, Lagoda, NEVallejo
# loci <- c(3:422)
inds <- c(1:33)
loci <- c(13:22)
# inds <- c(1:102) # all individuals
# loci <- c(3:4060) # ~ half the loci

small.df <- apo.df[inds, loci]

single.col.df <- small.df[, 0]
odd.cols <- seq(from = 1, to = ncol(small.df), by = 2)
locus <- 1
for (odd.col in odd.cols) {
  even.col <- odd.col + 1
  new.col <- paste(small.df[, odd.col], small.df[, even.col], sep = " ")
  col.name <- paste0("Locus", locus, sep = "")
  single.col.df[, col.name] <- new.col
  locus <- locus + 1
}

library("adegenet")
genid.2col <- df2genind(X = small.df,
                      ncode = 2,
                      ind.names = apo.df[inds, 1],
                      pop = apo.df[inds, 2],
                      ploidy = 2)

genid.1col <- df2genind(X = single.col.df,
                        sep = " ",
                        ncode = 2,
                        ind.names = apo.df[inds, 1],
                        pop = apo.df[inds, 2],
                        ploidy = 2)

# Make data frame with genotypes in a single column (not two columns)
# Start by replacing two-character genotypes with one-character genotypes
small.df.1char <- small.df
# replace 10 with 0, 11 with 1, etc...
replace.df <- data.frame(from = c(10, 11, 12, 13), to = c(0, 1, 2, 3))
for (i in 1:ncol(small.df.1char)) {
  for (j in 1:nrow(replace.df)) {
    matching.cells <- which(small.df.1char[, i] == replace.df[j, 1])
    small.df.1char[matching.cells, i] <- replace.df[j, 2]
  }
}
# Now combine values so genotype is in one, not two columns
# Also replace -9 with NA
small.df.1col <- small.df.1char[, 0]
odd.cols <- seq(from = 1, to = ncol(small.df.1char), by = 2)
locus <- 1
for (odd.col in odd.cols) {
  even.col <- odd.col + 1
  new.col <- paste0(small.df.1char[, odd.col], small.df.1char[, even.col])
  new.col <- gsub(pattern = "-9-9", replacement = NA, x = new.col)
  col.name <- paste0("Locus", locus, sep = "")
  small.df.1col[, col.name] <- new.col
  locus <- locus + 1
}
# Add pop id as first column
small.df.1col <- cbind(pop = apo.df[inds, 2], small.df.1col)

library("hierfstat")
p.fst.2col <- genet.dist(genid.2col, method = "WC84")
p.fst.1col <- genet.dist(genid.1col, method = "WC84")
p.fst.df <- genet.dist(dat = small.df.1col, diploid = TRUE, method = "WC84")

p.fst.2col
p.fst.1col


test.df <- data.frame(pop = c(1,1,1,2,2,2), one = c(00, 00, 00, 00, 00, 00), two = c(00, 10, 10, 00, 10, 11))

################################################################################
# Read in structure file as genid, then convert to hierfstat object
# genind2df
apo.str.genind <- read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str",
                                n.ind = 102,
                                n.loc = 4057,
                                onerowperind = TRUE,
                                col.lab = 1,
                                col.pop = 2,
                                col.others = 0,
                                row.marknames = 0,
                                NA.char = "-9",
                                sep = "\t")

# apo.str.df <- genind2df(x = apo.str.genind)

apo.str.hierfstat <- genind2hierfstat(dat = apo.str.genind)
start <- Sys.time()
p.fst <- genet.dist(apo.str.hierfstat, method="WC84") # ~10 minutes
end <- Sys.time()
# Still have really high values

################################################################################
# Copied & modified from Apo-hierfstat.R
# TODO: Compare the results of that code (looking at columns 3:433 [?]), to one 
# of the subsetting approaches above. 
# TODO: Also check about the issue of what looks like an incomplete genotype 
# reference...should be 3:432 - try that & compare to results in 
# output/Apo_pairwise_Fst.RData.  The numbers should not be wildly different
apo.df <- read.delim("data/Apodemia_0.9-noDockweiler-GNPSK.str", 
                     header = FALSE, 
                     na.strings = "-9", 
                     stringsAsFactors = FALSE)

# apo.df is:
# First column: sample name
# Second column: population (numeric)
# Remaining column: SNP data, with two columns per locus

###
# INCORRECT - READS EACH COLUMN AS A LOCUS
###

# See https://nescent.github.io/popgenInfo/DifferentiationSNP.html
library(adegenet)
library(hierfstat)
apo.genid <- df2genind(X = apo.df[,3:ncol(apo.df)], ploidy = 2, ind.names = apo.df[,1], pop = apo.df[,2], sep = "")
# Calculate pairwise Fst (takes ~21 minutes)
start <- Sys.time()
p.fst <- genet.dist(apo.genid, method = "WC84") # A 12 x 12 matrix (lower triangle only)
end <- Sys.time()
# TODO: Check to make sure # loci is correct! should be 4057 (half of 8116 - 2)
save(p.fst, file = "output/Apo-p-fst-na-strings.RData") 