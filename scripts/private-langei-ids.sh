#!/bin/bash

# Extract information about private alleles from VCF file
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-08-28

# Path information
ALLELEFREQS="../output/private-langei-allele-freqs.txt"
VCFFILE="../data/Apodemia_filteredVCF_0.9miss.recode.vcf"
OUTFILE="../output/private-langei-ids.txt"

# Start by grabbing the line from the VCF file that has column headers
VCFHEADER=$(grep -P 'CHROM\tPOS' $VCFFILE)
OUTHEAD=${VCFHEADER%QUAL*} # pull out stuff before QUAL
OUTHEAD=${OUTHEAD#'#'} # drop that leading #
echo -e $OUTHEAD > $OUTFILE

# Need to read in output/private-langei-allele-freqs.txt, line by line and
# pull out the value in the first column (skipping the first line)
readarray POSLINES < $ALLELEFREQS
LINES=0
for POSLINE in "${POSLINES[@]}";
do
  ((LINES++))
  if [ $LINES -ne 1 ];
  then
    # Pull out value in first column
    POS=$(echo $POSLINE | cut -d ' ' -f1)
    # Use that to grep the VCF file for the line we're interested in
    VCFLINE=$(grep -P 'un\t'${POS} $VCFFILE)
    # Drop everything after "PASS"
    KEEP=${VCFLINE%PASS*}
    # And kick out the last two characters (" .")
    KEEP2=${KEEP::-2}
    echo $KEEP2 >> $OUTFILE
  fi
done
