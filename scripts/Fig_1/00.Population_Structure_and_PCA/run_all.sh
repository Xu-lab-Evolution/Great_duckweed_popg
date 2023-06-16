data_dir="../../../data/" #diractory stores input VCF and temporary VCF
scripts_dir="./" #"scripts/Fig_1/00.Population_Structure_and_PCA/"
basic_set_SNP_vcf="SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.gz"
general_scripts_dir="../../other_scripts/"

#to filter basic SNP set with Hardy-Weinberg equilibrium test of P < 0.01 using vcftools
vcftools --gzvcf   \
	--hwe 0.01 \
	--not-chr pseudo0 \
	--remove-filtered-all \
	--recode --recode-INFO -all \
	--out $data_dir/SNP_vcf/SP_228.basic_set.snp.HWE001

#to prune SNP pairs with R-squre > 0.33 in a sliding window manner with 50 SNP window size and 5 SNP step size using Plinks
plink --vcf $data_dir/SNP_vcf/SP_228.basic_set.snp.HWE001.recode.vcf \
	--double-id --allow-extra-chr \
	--set-missing-var-ids @:# \
	--indep-pairwise 50 5 0.33 \
	--out $data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033

#to perform PCA analysis
plink $data_dir/SNP_vcf/SP_228.basic_set.snp.HWE001.recode.vcf \
	--double-id --allow-extra-chr --set-missing-var-ids @:# \
	--extract $data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.prune.in \
	--make-bed --pca \
	--out $data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.pca


#to perform population structure analysis using fastStructure, so please install it before use (https://github.com/rajanil/fastStructure)

for  i in $(seq 1 10); do python $general_scripts_dir/fastStructure/structure.py -K $i  --full  --seed=1234 --input=$data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.pca --output=$data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.simple  ; done

#to check best K value
python $general_scripts_dir/fastStructure/chooseK.py --input=$data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.simple > $data_dir/Fig_1/00.Population_Structure_and_PCA/chooseK.simple.res

#to visualize population structure
python  $general_scripts_dir/fastStructure/distruct.py -K 4 --title=K_4  --input=$data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.simple --output=$data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.simple.distruct.4.svg

#change color in SVG file
cp $data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.simple.distruct.4.svg $data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.simple.distruct.4.svg.fix.svg

fix_svg="$data_dir/Fig_1/00.Population_Structure_and_PCA/SP_228.hwe001.plk.ldRsqure0.33.no-pseudo0.simple.distruct.4.svg.fix.svg"
sed -i 's/\#b21212/\#FF0000/g'  $fix_svg
sed -i 's/\#62b212/\#39B54A/g'  $fix_svg
sed -i 's/\#12b2b2/\#F7931E/g'  $fix_svg
sed -i 's/\#6212b2/\#0000FF/g'  $fix_svg

#generating PCA_group_info

#PCA visualization using R, output plots will be generated in the "./output" dir
cd $scripts_dir/
Rscript SNP_PCA_plotting.R

#copy population structre plot (SVG) into output directory
cp $fix_svg ./output
