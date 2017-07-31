# Analyses for _Apodemia mormo_ project
## Last update: 2017-07-31

## Three folders included in the repo:
* data: the source data for analyses
* docs: documents, including outline, tasks, and early manuscript
* functions: R functions
* scripts: bash, python, and R scripts for various analyses
* tests: scripts primarily for testing functions under development

_**Note:**_ Before any R scripts are run, prepare the workspace:
1. create a directory called `output`
2. run the `scripts/prepare-data.sh` script to format data for all downstream
analyses

## data:
* Apo_localities.txt : a tab-delimited text file of localities with four 
columns: pop.name, pop.number, latitude, longitude
* Apodemia_0.9-noDockweiler-GNPSK.str : STRUCTURE-formatted (one-line) file 
excluding specimens from Dockweiler and Grasslands National Park populations 
(same data as in Apodemia_0.9miss.str sans specimens from those two populations)
* Apodemia_0.9miss.str : STRUCTURE-formatted (one-line) file, from JD including 
loci with up to 10% missing data (4057 loci)
* Apodemia_filteredVCF_0.9miss.recode.vcf : a VCF file that was precursor to JD 
STRUCTURE runs (10% missing data allowed) (4057 loci)
* inds-all.txt : text file with list of all individual identifiers (e.g. 
ArroyoBayo_26, langei_JO1671)
* pop-prefixes.txt : text file with substrings corresponding to population 
identifying portion of individual identifiers (e.g. ArroyoBayo, langei)

***

## docs:
* apo-ms-outline.md : early outline of manuscript
* apo-ms-to-doc.sh : bash script for converting markdown `apo-ms.md` to a 
Word document
* apo-ms.md : early draft of manuscript
* apo-pandoc-ref.docx : a Word document to serve as a template for pandoc (called 
in `apo-ms-to-doc.sh`) when converting `apo-ms.md` to a Word document
* apo-task-list.md : a task list

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

## tests
* rarefaction-test-01.R
* rarefaction-test-02.R
* rarefaction-test-03.R
* rarefaction-test-04.R
* rarefaction-test-05.R
* rarefaction-test-06.R
* rarefaction-test-07.R
* rarefaction-test-08.R
* rarefaction-test-09.R
* rarefaction-test-10.R
* rarefaction-test-11.R
* rarefaction-test-12.R
* rarefaction-test-13.R
* rarefaction-test-14.R
