# Partial code and data for Dupuis et al. 2018
[https://doi.org/10.1007/s10592-018-1081-8](https://doi.org/10.1007/s10592-018-1081-8)

Julian R. Dupuis, Jeffrey C. Oliver, Bryan M.T. Brunet, Travis Longcore, Jana J. Johnson, and Felix A.H. Sperling. 2018. Genomic data indicate ubiquitous evolutionary distinctiveness among populations of California metalmark butterflies. _Conservation Genetics_.

## Three folders included in the repo:
* data: the source data for analyses
* functions: R functions
* scripts: bash, python, and R scripts for various analyses

_**Note:**_ Before any R scripts are run, prepare the workspace:
1. create a directory called `output`
2. run the `scripts/prepare-data.sh` script to format data for all downstream
analyses

## data:
* Apo_localities.txt : a tab-delimited text file of localities with four
columns: pop.name, pop.number, latitude, longitude
* Apodemia_0.9-noDockweiler-GNPSK.str : STRUCTURE-formatted (one-line) file
excluding specimens from Dockweiler and Grasslands National Park populations;
same data as in Apodemia_0.9miss.str sans specimens from those two populations,
produced via `grep -v "Dockweiler" Apodemia_0.9miss.str | grep -v "NPSK" > Apodemia_0.9-noDockweiler-GNPSK.str`
* Apodemia_0.9-noGNPSK.str : STRUCTURE-formatted (one-line) file excluding specimens
from Grasslands National Park population; same data as in Apodemia_0.9miss.str sans
specimens from Grasslands NP, produced via `grep -v grep -v "NPSK" > Apodemia_0.9-noGNPSK.str`
* Apodemia_0.9miss.str : STRUCTURE-formatted (one-line) file, from JD including
loci with up to 10% missing data (4057 loci)
* Apodemia_filteredVCF_0.9miss.recode.vcf : a VCF file that was precursor to JD
STRUCTURE runs (10% missing data allowed) (4057 loci)
* inds-all.txt : text file with list of all individual identifiers (e.g.
ArroyoBayo_26, langei_JO1671)
* pop-prefixes.txt : text file with substrings corresponding to population
identifying portion of individual identifiers (e.g. ArroyoBayo, langei)

***

## functions:
* apodemia-functions.R : functions written to make the code in scripts a bit
more readable
* rarefaction-functions.R : code for running private allele and allelic richness
rarefaction to control for sampling efforts
***

## scripts:
* Apo-fastStructure-template.sh : executes multiple runs of STRUCTURE
* ibd-plot.R : scatterplots of genetic distance vs. geographic distance
* ibd.R : tests for isolation by distance
* mapping.R : map of populations
* pairwise-fst-boxplot.R : boxplot of pairwise Fst values for each population
* pairwise-fst-vcf.sh : uses vcftools to calculate pairwise Fst values; used
primarily as a reality check on values from the hierfstat package genet.dist
function (weighted Fst estimate from vcftools returns results identical to
those produced by genet.dist(method = "WC85"))
* prepare-data.R : various data preparation for downstream analyses
* rarefaction-richness-private-permutation.R : Estimates allelic richeness and
private allele counts for each population, taking sampling effort into account.
For more details on how data are rarefied, see comments in `functions/rarefaction-functions.R`
* str_to_two_line.py : convert structure-formatted file from single-line to
two-line format; primarily to create files for use with fastStructure