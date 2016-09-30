#!/bin/bash

# For converting a nexus file to phylip
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-03-02

# Install Biopython
# sudo apt-get install python-biopython
# Install pip for managing python packages
# sudo apt-get install python-pip
# Then install seqmagick
# sudo pip install seqmagick
# Convert to 'relaxed' Phylip format
seqmagick convert ../Data/Apodemia_0.5miss.nex ../Data/Apodemia-2016-05-12.phyx

