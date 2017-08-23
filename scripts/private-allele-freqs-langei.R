# Script for identifying langei private alleles
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-07-28

################################################################################
# SUMMARY
# *

################################################################################
# SETUP
# Identify data files
langei.file <- "output/allele-freqs/allele-freqs-langei.frq"
non.langei.file <- "output/allele-freqs/allele-freqs-non-langei.frq"

################################################################################
# DATA PREP
# Read in data, update column names, and combine frequencies into one data frame

# Read in data
langei.freqs <- read.delim(file = langei.file, header = TRUE)
non.langei.freqs <- read.delim(file = non.langei.file, header = TRUE)

# Extract allele frequencies & marker id
combined.freqs <- langei.freqs[, c("POS", "REF", "ALT")]

# Update column names
colnames(combined.freqs)[which(colnames(combined.freqs) == "POS")] <- "pos"
colnames(combined.freqs)[which(colnames(combined.freqs) == "REF")] <- "langei.ref"
colnames(combined.freqs)[which(colnames(combined.freqs) == "ALT")] <- "langei.alt"

# Calculate number of individuals sampled (N_CHR is number of chromosomes, so 
# N_CHR = 22 means 11 individuals were sampled)
combined.freqs$n.langei <- langei.freqs$N_CHR / 2

# Merge data frames based on marker id (pos)
colnames(non.langei.freqs)[which(colnames(non.langei.freqs) == "REF")] <- "non.langei.ref"
colnames(non.langei.freqs)[which(colnames(non.langei.freqs) == "ALT")] <- "non.langei.alt"
combined.freqs <- merge(x = combined.freqs, 
                        y = non.langei.freqs[, c("POS", "non.langei.ref", "non.langei.alt")],
                        by.x = "pos",
                        by.y = "POS")

################################################################################
# ANALYZE
# Identify private alleles, subset data accordingly, and write to file

# Identify those alleles fixed for the reference allele in non-langei 
# populations
fixed.ref.non.langei <- combined.freqs[which(combined.freqs$non.langei.alt == 0), ]

# TODO: what about markers with alternate allele fixed in non.langei 
# populations? Doesn't happen with these data, but could miss private allele if 
# the langei REF allele is only found in langei

write.table(x = fixed.ref.non.langei, 
            file = "output/private-langei-allele-freqs.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
