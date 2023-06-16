#module load GCC/8.2.0-2.31.1 BEDTools/2.28.0
#ml bio/SAMtools/1.14-GCC-11.2.0 #@ mogon cluster, change it when naccesary
samtools faidx ../../../data//Spolyrhiza_ref/SP_combined.fasta
awk '{OFS="\t"; print $1,$2}' ../../../data//Spolyrhiza_ref/SP_combined.fasta.fai > ../../../data//Spolyrhiza_ref/SP_combined.fasta.genome.txt
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$1."_"$2,$4"/"$5}}' ../../../data//SNP_vcf/SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf > ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.bed

#ml bio/BEDTools/2.29.2-iccifort-2020.1.217
bedtools slop -i ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.bed  -g ../../../data//Spolyrhiza_ref/SP_combined.fasta.genome.txt -l 150 -r 150  > ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.flanking150.bed

awk '{if($3-$2==301){print}}' ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.flanking150.bed | cat > ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.flanking150.fixed.bed

bedtools getfasta  -fi ../../../data//Spolyrhiza_ref/SP_combined.fasta  -name -bed ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.flanking150.fixed.bed > ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.flanking150.fa
