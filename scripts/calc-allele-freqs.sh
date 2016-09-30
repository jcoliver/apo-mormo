# For each population, calculate allele frequencies for that population and all 
# remaning populations. e.g. calculate langei frequences and all remaining non-
# langei populations

# See original implementation in calc-allele-freqs-langei.sh

while read POP; do
	echo "========== = = $POP = = =========="

	# text file that includes only those specimens from POP
	grep -e "$POP" inds-all.txt > inds-$POP.tmp

	# Pull out allele frequencies for POP
	vcftools --vcf Apodemia-all-pops.vcf --keep inds-$POP.tmp --recode --stdout | vcftools --vcf - --freq2 --out allele-freqs-$POP
	rm allele-freqs-$POP.log

	# Replace {FREQ} header with more informative REF/ALT
	sed -i 's/{FREQ}/REF\tALT/g' allele-freqs-$POP.frq
	# Convert -nan to R standard "NA"
	sed -i 's/-nan/NA/g' allele-freqs-$POP.frq

	# Move to directory for R analyses
	mv allele-freqs-$POP.frq ../Apo-R/data/allele-freqs-$POP.frq

	# Calculate allele frequences for all specimens not in POP
	vcftools --vcf Apodemia-all-pops.vcf --remove inds-$POP.tmp --recode --stdout | vcftools --vcf - --freq2 --out allele-freqs-non-$POP
	rm allele-freqs-non-$POP.log

	# Replace {FREQ} header with more informative REF/ALT
	sed -i 's/{FREQ}/REF\tALT/g' allele-freqs-non-$POP.frq
	# Convert -nan to R standard "NA"
	sed -i 's/-nan/NA/g' allele-freqs-non-$POP.frq

	# Move to directory for R analyses
	mv allele-freqs-non-$POP.frq ../Apo-R/data/allele-freqs-non-$POP.frq

	rm inds-$POP.tmp
done < pop-prefixes.txt

