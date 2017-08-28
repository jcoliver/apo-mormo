#!/bin/bash

# Extract information about private alleles from VCF file
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-08-28

# Need to read in output/private-langei-allele-freqs.txt, line by line and
# pull out the value in the first column

# For each column, want to grep on that value for line in VCF file
# where mynum is the value in the pos column
grep -P 'un\t'${mynum} data/Apodemia_filteredVCF_0.9miss.recode.vcf 

# Pull out some number characters(?) from that output

# When using grep, need to use the -P flag to get it to recognize tab \t