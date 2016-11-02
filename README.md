# Analyses for _Apodemia mormo_ project
## Last update: 2016-11-01

## Three folders included in the repo:
1. data: the source data for analyses
2. functions: R functions
3. scripts: bash, python, and R scripts for various analyses

## 1. data:
* Apo_localities.txt : a tab-delimited text file of localities with four 
columns: pop.name, pop.number, latitude, longitude
* Apo_0.9-compare-str-df.str : ????
* Apodemia_0.9-for-R.str : STRUCTURE-formatted file for use in R work
* Apodemia_0.9-full-test.str : ????
* Apodemia_0.9-medium-test.str : ????
* Apodemia_0.9-noDockweiler-GNPSK.str : STRUCTURE-formatted (one-line) file 
excluding specimens from Dockweiler and Grasslands National Park populations 
(same data as in Apodemia_0.9miss.str sans specimens from those two populations)
* Apodemia_0.9-small-test.str : ????
* Apodemia_0.9miss.str : STRUCTURE-formatted (one-line) file, from JD including 
loci with up to 10% missing data (4057 loci)
* Apodemia_filteredVCF_0.9miss.recode.vcf : a VCF file that was precursor to JD 
STRUCTURE runs (10% missing data allowed) (4057 loci)
* filtered_Apo.str : ????
* filtered_noGNPSK.str : ????
* inds-all.txt : text file with list of all individual identifiers (e.g. 
ArroyoBayo_26, langei_JO1671)
* pop-prefixes.txt : text file with substrings corresponding to population 
identifying portion of individual identifiers (e.g. ArroyoBayo, langei)

***

## 2. functions:
* apodemia-functions.R : functions written to make the code in scripts a bit 
more readable

***

## 3. scripts:
* allelic-richness.R : INCOMPLETE calculates the allelic richness (# of alleles) 
in each population
* Apo-fastStructure-template.sh : executes multiple runs of STRUCTURE
* apo-nexus-to-phylip.sh : converts nexus to phylip format; relies on seqmagick
* apo-raxml-commands.sh : RAxML bootstrapping
* calc-allele-freqs.sh : calculates allele frequencies for each population
* create-langei-fixed-private-vcf.sh : creates a VCF file containing only those
loci which are private to and fixed in langei population
* dapc.R : INCOMPLETE discriminant analysis of principal components
* ibd-plot.R : scatterplots of genetic distance vs. geographic distance
* ibd.R : tests for isolation by distance
* mapping.R : map of populations
* pairwise-fst-boxplot.R : boxplot of pairwise Fst values for each population
* pairwise-fst-vcf.sh : uses vcftools to calculate pairwise Fst values; used 
primarily as a reality check on values from the hierfstat package genet.dist 
function (weighted Fst estimate from vcftools returns results identical to 
those produced by genet.dist(method = "WC85"))
* pca.R : principal components analysis plots
* prepare-data.R : various data preparation for downstream analyses
* private-allele-counts.R : INCOMPLETE counts the number of private alleles for 
each population
* private-allele-freqs-langei.R : INCOMPLETE creates table of allele frequencies
for all loci that have variants only found in the langei population
