#!/bin/bash

# Create VCF file with only those alleles unique to langei that are fixed
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-30

cd ~/Documents/Other\ Work/Apodemia/apo-repo/scripts
VCFFILE="../data/Apodemia_filteredVCF_0.9miss.recode.vcf"
OUTFILE="../output/langei-private-fixed.vcf"
HEADERLINE="$(head $VCFFILE | tail -n1)"
ALLELES="$(grep -e $'\t''128466'$'\t' -e $'\t''200791'$'\t' -e \
$'\t''845685'$'\t' -e $'\t''1320698'$'\t' $VCFFILE)"

echo "${HEADERLINE}" > $OUTFILE
echo "${ALLELES}" >> $OUTFILE
