#!/bin/bash

# Create file with just VCF header line
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-30

VCFFILE="../data/Apodemia_filteredVCF_0.9miss.recode.vcf"
OUTFILE="../output/Apodemia-vcf-header.txt"

head $VCFFILE | tail -n1 > $OUTFILE
