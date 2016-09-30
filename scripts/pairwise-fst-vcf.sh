# Do all pairwise-Fst comparisons among populations
# Execute from ~/Documents/Other Work/Apodemia/Code
# Note this pulls out the weighted Fst estimate, which appears to 
# be the same value calculated by the genet.dist(method = "WC85")
# function in the hierfstat package for R

VCFFILE="../Data/Apodemia_filteredVCF_0.9miss.recode.vcf"
RESULTFILE="../Output/pairwise-fst-vcf.txt"
POPPREFIX="../Data/pop-prefixes.txt"
ALLINDS="../Data/inds-all.txt"
FSTTEMP="fst.tmp"

#POP1="HullMtn"
#POP2="langei"
#cat $POPPREFIX > pop1-prefixes.txt
OUTERPOPS="$(cat $POPPREFIX)"
SKIPPOP="$(head -n1 $POPPREFIX)"
grep -v "$SKIPPOP" $POPPREFIX > innerpops.tmp
echo -e "POP1\tPOP2\tFst" > $RESULTFILE

for POP1 in $OUTERPOPS; do
	POP1FILE="inds-$POP1.tmp"
	grep -e "$POP1" $ALLINDS > $POP1FILE
	INNERPOPS="$(cat innerpops.tmp)"
	FIRSTPOP2="$(head -n1 innerpops.tmp)"
	for POP2 in $INNERPOPS; do
		echo "===   $POP1 vs. $POP2   ==="
		POP2FILE="inds-$POP2.tmp"
		grep -e "$POP2" $ALLINDS > $POP2FILE
		# output by vcftools sent to stderr instead of stdout, hence 2> instead of >
		vcftools --vcf $VCFFILE --weir-fst-pop $POP1FILE --weir-fst-pop $POP2FILE 2> $FSTTEMP
		#FSTLINE="$(grep -e "Weir and Cockerham mean Fst estimate:" "$FSTTEMP")"
		FSTLINE="$(grep -e "Weir and Cockerham weighted Fst estimate:" "$FSTTEMP")"
		FST="$(echo ${FSTLINE##*: })"
		echo -e "$POP1\t$POP2\t$FST" >> $RESULTFILE
		rm $POP2FILE $FSTTEMP
	done
	rm $POP1FILE 
	grep -v "$FIRSTPOP2" innerpops.tmp > new-innerpops.tmp
	mv new-innerpops.tmp innerpops.tmp
done

rm innerpops.tmp
rm out.weir.fst

# Development
# POP1="HullMtn"
# POP2="langei"
# POP1FILE="inds-$POP1.tmp"
# POP2FILE="inds-$POP2.tmp"
# ALLINDS="../Data/inds-all.txt"
# grep -e "$POP1" $ALLINDS > $POP1FILE
# grep -e "$POP2" $ALLINDS > $POP2FILE
# vcftools --vcf ../Data/Apodemia_filteredVCF_0.9miss.recode.vcf --weir-fst-pop $POP1FILE --weir-fst-pop $POP2FILE
# rm out.log out.weir.fst $POP1FILE $POP2FILE
