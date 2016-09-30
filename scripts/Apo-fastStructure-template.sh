#!/bin/bash

# INCOMPLETE: suspended development as J. Dupuis is handling STRUCTURE runs

# Running fastStructure on Apodemia data
# using file created by str_to_two_line.py script

#TODO: Need to parameterize this:
# $fspath (path to fastSTRUCTURE)
# $infile
# $outfile (will be the input for chooseK.py)
# $rngseed
# $popfile
# output filenames for graphics should just be:
#   $outfile ".distruct." K ".svg (e.g. Apo/ApoTest05.distruct.2.svg)

## Reporting function if parameters are not correct
usage() {
	echo "Usage: $0 -i <input file> -o <output file> -k <max K> -p <popfile> [-s <seed>]"
	exit 1
}


# in ~/Executables/fastStructure:
python structure.py -K 1 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 2 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 3 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 4 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 5 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 6 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 7 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 8 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 9 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 10 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 11 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 12 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 13 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 14 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str
python structure.py -K 15 --input=Apo/Apo_fastStr3 --output=Apo/ApoTest05 --full --seed=201511251522 --format=str

# Compare models (it only looks at Marginal Likelihoods, not the delta statistic):
python chooseK.py --input=Apo/ApoTest05

# Model complexity that maximizes marginal likelihood = 2
# Model components used to explain structure in data = 5

# Graphics
# For --input, use same value as used in structure.py call
# Need to create a popfile, see create_popfile.py script
python distruct.py -K 2 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.2.svg --popfile=Apo/Apo.popfile --title="K = 2"
python distruct.py -K 3 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.3.svg --popfile=Apo/Apo.popfile --title="K = 3"
python distruct.py -K 4 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.4.svg --popfile=Apo/Apo.popfile --title="K = 4"
python distruct.py -K 5 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.5.svg --popfile=Apo/Apo.popfile --title="K = 5"
python distruct.py -K 6 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.6.svg --popfile=Apo/Apo.popfile --title="K = 6"
python distruct.py -K 7 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.7.svg --popfile=Apo/Apo.popfile --title="K = 7"
python distruct.py -K 8 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.8.svg --popfile=Apo/Apo.popfile --title="K = 8"
python distruct.py -K 9 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.9.svg --popfile=Apo/Apo.popfile --title="K = 9"
python distruct.py -K 10 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.10.svg --popfile=Apo/Apo.popfile --title="K = 10"
python distruct.py -K 11 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.11.svg --popfile=Apo/Apo.popfile --title="K = 11"
python distruct.py -K 12 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.12.svg --popfile=Apo/Apo.popfile --title="K = 12"
python distruct.py -K 13 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.13.svg --popfile=Apo/Apo.popfile --title="K = 13"
python distruct.py -K 14 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.14.svg --popfile=Apo/Apo.popfile --title="K = 14"
python distruct.py -K 15 --input=Apo/ApoTest05 --output=Apo/ApoTest05.distruct.15.svg --popfile=Apo/Apo.popfile --title="K = 15"

