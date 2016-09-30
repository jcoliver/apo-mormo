# Calculating Fst to compare with R results
# Executed from ~/Documents/Other Work/Apodemia/Code

# Current data
vcftools --vcf ../Data/Apodemia_filteredVCF_0.9miss.recode.vcf --weir-fst-pop ../Data/inds-langei.txt --weir-fst-pop ../Data/inds-HullMtn.txt
vcftools --vcf ../Data/Apodemia_filteredVCF_0.9miss.recode.vcf --weir-fst-pop ../Data/inds-langei.txt --weir-fst-pop ../Data/inds-Ladoga.txt
vcftools --vcf ../Data/Apodemia_filteredVCF_0.9miss.recode.vcf --weir-fst-pop ../Data/inds-langei.txt --weir-fst-pop ../Data/inds-NEVallejo.txt
vcftools --vcf ../Data/Apodemia_filteredVCF_0.9miss.recode.vcf --weir-fst-pop ../Data/inds-ArroyoBayo.txt --weir-fst-pop ../Data/inds-DelPeurtoCanyon.txt

# Previous, smaller data set
vcftools --vcf ../Data/filtered_loci_95missing.vcf --weir-fst-pop ../Data/inds-langei.txt --weir-fst-pop ../Data/inds-HullMtn.txt
vcftools --vcf ../Data/filtered_loci_95missing.vcf --weir-fst-pop ../Data/inds-langei.txt --weir-fst-pop ../Data/inds-Ladoga.txt
vcftools --vcf ../Data/filtered_loci_95missing.vcf --weir-fst-pop ../Data/inds-langei.txt --weir-fst-pop ../Data/inds-NEVallejo.txt
vcftools --vcf ../Data/filtered_loci_95missing.vcf --weir-fst-pop ../Data/inds-ArroyoBayo.txt --weir-fst-pop ../Data/inds-DelPeurtoCanyon.txt

rm out.weir.fst





