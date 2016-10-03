# Script for identifying langei-unique SNPs
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-07-28

rm(list = ls())
# Read in vcftools-produced frequency files
langei.freqs <- read.delim(file = "data/allele-freqs-langei.frq", header = TRUE)
rownames(langei.freqs) <- paste0("pos", langei.freqs$POS)
non.langei.freqs <- read.delim(file = "data/allele-freqs-non-langei.frq", header = TRUE)

combined.freqs <- langei.freqs[, c("REF", "ALT")]
colnames(combined.freqs)[which(colnames(combined.freqs) == "REF")] <- "langei.ref"
colnames(combined.freqs)[which(colnames(combined.freqs) == "ALT")] <- "langei.alt"
combined.freqs$n.langei <- langei.freqs$N_CHR / 2

combined.freqs$non.langei.ref <- non.langei.freqs$REF
combined.freqs$non.langei.alt <- non.langei.freqs$ALT

no.alt.in.non.langei <- combined.freqs[which(combined.freqs$non.langei.alt == 0), ]
fixed.alt.langei <- no.alt.in.non.langei[which(no.alt.in.non.langei$langei.ref == 0), ]

# Saving some data for identifying genes
no.alt.in.non.langei$pos <- rownames(no.alt.in.non.langei)
no.alt.in.non.langei$pos <- gsub("pos", "", no.alt.in.non.langei$pos)
write.table(x = no.alt.in.non.langei, file = "output/langei-unique-snps-info.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
