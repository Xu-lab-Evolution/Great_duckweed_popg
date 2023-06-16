data_dir="../../data/" #diractory stores input VCF and temporary VCF
scripts_dir="./" #"scripts/Fig_1/00.Population_Structure_and_PCA/"
basic_set_SNP_vcf="SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.gz"
general_scripts_dir="../other_scripts/"
pi_dir="$data_dir/Fig_2/00.pi_and_piN_piS/"
het_dir="$data_dir/Fig_2/01.genome_wide_het/"

#################################################00.caculate genome-wide pi and piN/piS
#to filter out SNPs that have missing data
vcftools --vcf $data_dir/SNP_vcf/SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf --max-missing 1.0  --remove-filtered-all --recode --recode-INFO-all --out $data_dir/SNP_vcf/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf

#to filter SNPs that overlapped with SV
#ml bio/BEDTools/2.29.2-iccifort-2020.1.217 #@load bedtools module on the mogon cluster
bedtools subtract -header -a $data_dir/SNP_vcf/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.recode.vcf  -b $data_dir/SV_vcf/SV.bed > $pi_dir/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.rm_SV.vcf
#ml bio/BCFtools/1.14-GCC-11.2.0 #@load bedtools module on the mogon cluster
bgzip $pi_dir/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.rm_SV.vcf && tabix -p vcf $pi_dir/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.rm_SV.vcf.gz

#extract each population from VCF
#ml bio/BCFtools/1.14-GCC-11.2.0 #@load bedtools module on the mogon cluster
for i in AME ASIA EUR IND; do mkdir -p $pi_dir/$i ; grep $i $data_dir/clonal_family/family_info.VCF_ID.txt | awk '{print $5}' > $pi_dir/$i.list ; bcftools view -S $pi_dir/$i.list $pi_dir/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.rm_SV.vcf.gz -Oz -o $pi_dir/$i/$i.vcf.gz & done

#split each population's VCF based on each chromosome
mkdir -p $pi_dir/ALLPOP ; cp $pi_dir/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.rm_SV.vcf.gz $pi_dir/ALLPOP/ALLPOP.vcf.gz ; cp $pi_dir/SP_228.basic_set.snp.recode.rm_cluster3-10.0missing.vcf.rm_SV.vcf.gz.tbi $pi_dir/ALLPOP/ALLPOP.vcf.gz.tbi  

for i in AME ASIA EUR IND ALLPOP; do $general_scripts_dir/split_VCF_based_on_each_chr.sh $pi_dir/$i/$i.vcf.gz; for chr in `cat $pi_dir/$i/chromosomes.txt`; do mkdir -p $pi_dir/$i/$chr ; mv $pi_dir/$i/$chr.vcf $pi_dir/$i/$chr ; done ; done

#first correct the allele frequency information (AF tag) from each VCF file, then perform SNPgenie calculation
if [ -f snpgenie.all_run.list.sh ]; then rm snpgenie.all_run.list.sh ; fi
#for i in AME ASIA EUR IND ALLPOP; do chr=(ChrS{01..20}); for chr in ${chr[@]}; do perl get_snpgenie_cmd.pl $PWD/$pi_dir/$i/$chr/$chr.vcf $PWD/$data_dir/Spolyrhiza_ref/split_by_chr $PWD/$data_dir/Spolyrhiza_annotation_SPGA2022/split_chr $PWD/$general_scripts_dir/SNPgenie/ mogon ; done ; done #run on mogon
for i in AME ASIA EUR IND ALLPOP; do chr=(ChrS{01..20}); for chr in ${chr[@]}; do perl get_snpgenie_cmd.pl $PWD/$pi_dir/$i/$chr/$chr.vcf $PWD/$data_dir/Spolyrhiza_ref/split_by_chr $PWD/$data_dir/Spolyrhiza_annotation_SPGA2022/split_chr $PWD/$general_scripts_dir/SNPgenie/ ; done ; done

#collect result from population_summary.txt (SNPgenie's result)
find $pi_dir  -name "population_summary.txt" | sort > $pi_dir/SNPgenie_res.list

perl dealing_Res.pl $pi_dir/SNPgenie_res.list

mkdir -p ./output ; cp $pi_dir/SNPGenie_res_dealing.res ./output 

#################################################01.caculate the heterozygosity
vcftools --vcf $data_dir/SNP_vcf/SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf --het  --out $het_dir/output_het.MAF0.05  --maf 0.05
Rscript het_box_plot.R $data_dir/clonal_family/family_info.pop_full_name.txt  $het_dir/output_het.MAF0.05.het

cp $het_dir/MAF0_05.het_R.pdf ./output 
cp $het_dir/het_df.MAF_0_05.tsv ./output 

#################################################02.caculate the genome-wide recombination rate




#################################################summarize and visualization
#ml lang/R/4.2.0-foss-2021b #@mogon 
Rscript 4pop_2SV.piechart.R
Rscript 5_parameters.barplot.R

