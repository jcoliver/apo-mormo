# Script for identifying population-unique SNPs
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-07-29

rm(list = ls())
pops <- scan(file = "data/pop-prefixes.txt", what = character())

results.df <- data.frame(pop = pops, unique.snps = NA, fixed.unique.snps = NA)

for (one.pop in pops) {
  focal.freqs <- read.delim(file = paste0("output/allele-freqs/allele-freqs-", one.pop, ".frq"), header = TRUE)
  other.freqs <- read.delim(file = paste0("output/allele-freqs/allele-freqs-non-", one.pop, ".frq"), header = TRUE)
  
  combined.freqs <- focal.freqs[, c("POS", "CHROM")]
  combined.freqs$focal.ref <- focal.freqs$REF
  combined.freqs$focal.alt <- focal.freqs$ALT
  combined.freqs$other.ref <- other.freqs$REF
  combined.freqs$other.alt <- other.freqs$ALT
  
  unique.focal.alleles <- combined.freqs[which(combined.freqs$other.alt == 0), ]

  if (sum(combined.freqs$other.ref == 0) > 0) {
    cat("ALT allele fixed in non-", one.pop, " populations.\n", sep = "")
  }
  
  fixed.unique.focal.alleles <- unique.focal.alleles[which(unique.focal.alleles$focal.ref == 0), ]
  
  results.row <- which(results.df$pop == one.pop)
  results.df$unique.snps[results.row] <- nrow(unique.focal.alleles)
  results.df$fixed.unique.snps[results.row] <- nrow(fixed.unique.focal.alleles)

}

date.for.filenames <- format(Sys.Date(), "%Y-%m-%d")
write.table(x = results.df, file = paste0("output/unique-snps-", date.for.filenames, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
