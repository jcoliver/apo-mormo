#!/usr/bin/env python2.7

# Converts a single-line structure file into 
# two-line format required by fastStructure
# Input file format:
# Sample_01 1   10  10  13  11  12  12
# Sample_02 1   12  12  -9  -9  12  12
# Sample_03 2   12  12  13  13  10  12
# This input has PAIRS of alleles in successive
# columns, with bases corresponding to a two-
# digit integer (A=10,T=11,G=12,C=13), so the 
# genotypes for the example above would be:
# Sample_01 1   AA  CT  GG
# Sample_02 1   GG  -9  GG
# Sample_03 2   GG  CC  AG

# Output file format: (breaking each sample 
# into two lines and adding four empty columns
# [spaces] at the beginning of each line)
# bases now represented by a single-digit 
# integer (A=0,T=1,G=2,C=3)
#     Sample_01 1 0 3 2
#     Sample_01 1 0 2 2
#     Sample_02 1 2 -9 2
#     Sample_02 1 2 -9 2
#     Sample_03 2 2 3 0
#     Sample_03 2 2 3 2

import io
import sys
import os.path

class Sample:
    def __init__(self, sample_name, population, genotypes):
        self.name = sample_name
        self.population = population
        self.genotypes = genotypes

# Check for str file
args = sys.argv
if len(args) < 3:
    print "Usage str_to_two_line.py <str file> <output filename>"
    # str file should be something like filtered_Apo.str
    # output filename should be something like Apo_fastStr3.str
    quit()

# Make sure file exists
str_filename = args[1]
if not(os.path.isfile(str_filename)):
    print "str file \"" + str_filename + "\" does not exist"
    quit()

# Make sure output file does not already exist
out_filename = args[2]
if os.path.isfile(out_filename):
    print "output filename \"" + out_filename + "\" already in use"
    quit()

# Open the pipe to file and read each line, splitting
# each genoptype into two alleles...
# Samples [dict]
#   loci [list]
#       alleles [list]

samples = list()
str_fh = io.open(str_filename, "rb")
hit_warning_limit = False
warned_count = 0

allele_conv = dict()
allele_conv["10"] = 0 # A
allele_conv["11"] = 1 # T
allele_conv["12"] = 2 # G
allele_conv["13"] = 3 # C
allele_conv["-9"] = -9 # Missing

for line in str_fh:
    linestripped = line.strip()
    if len(linestripped) > 0:
        linesplit = linestripped.split()
        sample_name = linesplit[0]
        population = linesplit[1]
        # A list of 2-element lists
        genotypes = list()
#        for i in range(2, len(linesplit)):
        for i in range(2, len(linesplit), 2):
            allele_one = linesplit[i]
            allele_two = linesplit[i + 1]
            if not(allele_one in allele_conv.keys()):
                print "Unrecognized allele one in " + str(sample_name) + ": " + str(allele_one)
                quit()
            if not(allele_two in allele_conv.keys()):
                print "Unrecognized allele two in " + str(sample_name) + ": " + str(allele_two)
                quit()
            alleles = [allele_conv[allele_one], allele_conv[allele_two]]
#            genotype = linesplit[i]
#            alleles = [-9,-9]
#            if genotype != "-9":
#                alleles = [genotype[0], genotype[1]]
            genotypes.append(alleles)
        if len(genotypes) > 0:
            new_sample = Sample(sample_name, population, genotypes)
            samples.append(new_sample)
        else:
            if not(hit_warning_limit):
                print "Zero data for sample " + str(sample_name)
                warned_count = warned_count + 1
                if (warned_count >= 5):
                    hit_warning_limit = True
    # End conditional for non-empty line
# End looping over all lines

# use this to start a line
four_spaces = "    "
sep = " "

outfile_fh = io.open(out_filename, "wb")
# loop over the dictionary and write to the file
first_line = True
for sample in samples:
    if first_line:
        one_line = ""
        first_line = False
    else:
        one_line = "\n"

    line_start = four_spaces + str(sample.name) + sep + str(sample.population)
    one_line = one_line + line_start
    two_line = "\n" + line_start
    for genotype in sample.genotypes:
        one_line = one_line + sep + str(genotype[0])
        two_line = two_line + sep + str(genotype[1])

    outfile_fh.write(one_line)
    outfile_fh.write(two_line)
# End looping over samples_to_loci dict
outfile_fh.close()

