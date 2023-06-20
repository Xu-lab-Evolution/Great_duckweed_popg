vcftools=/home/pablo/Programs/vcftools/src/cpp/vcftools
vcfAME=/home/pablo/Documents/Duckweed/original_vcfs/AMEbcftools.vcf
vcfAll=/home/pablo/Documents/Duckweed/original_vcfs/all_chroms_renamed_ordered_selected.recode.vcf

posAME=/home/pablo/Documents/Duckweed/original_vcfs/positionsAME_polymorphic.txt
posAME_chr_01=/home/pablo/Documents/Duckweed/original_vcfs/positionsAME_polymorphic_chr_01.txt

$vcftools --vcf $vcfAll --positions $posAME --keep ASIA.txt --counts2 --out countsASIA
$vcftools --vcf $vcfAll --positions $posAME --keep AME.txt --counts2 --out countsAME
$vcftools --vcf $vcfAll --positions $posAME --keep EUR.txt --counts2 --out countsEUR
$vcftools --vcf $vcfAll --positions $posAME --keep IND.txt --counts2 --out countsIND

#$vcftools --vcf $vcfAME --counts2 --out countsAME



