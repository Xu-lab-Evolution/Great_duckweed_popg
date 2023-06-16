data_dir="../../../data/" #diractory stores input VCF and temporary VCF
scripts_dir=`pwd` #"scripts/Fig_1/00.Population_Structure_and_PCA/"
general_scripts_dir="../../other_scripts/"
basic_set_SNP_vcf="SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.gz"
ld_dir="$data_dir/Fig_1/02.Linkage_disequilibrium/LD/"

mkdir -p $ld_dir
mkdir -p $ld_dir/split_pop;

#to split VCF into diferent populations, and also to keep one genotype from each clonal family
for i in "AME" "ASIA" "EUR" "IND"; do echo $i ; mkdir -p $ld_dir/split_pop/$i/ ; bcftools view -Ov -S $data_dir/clonal_family/$i.remove_clonality.txt -o $ld_dir/split_pop/$i/SP228.basic_set.$i.vcf $data_dir/SNP_vcf/$basic_set_SNP_vcf ; done

#to further split population VCF to each chromosome's VCF
#ml bio/BCFtools/1.14-GCC-11.2.0 #load module @ mogon
for i in "AME" "ASIA" "EUR" "IND"; do echo $i ;  $general_scripts_dir/split_VCF_based_on_each_chr.sh $ld_dir/split_pop/$i/SP228.basic_set.$i.vcf & done

#ml lib/zlib/1.2.12 #load module @ mogon
#to generate PopLDdecay running script (split based on each chromesome and per population)
#let's assume the PopLDdecay has been properly installed and the executable was put in the default path

if [ -f "PopLDdecay_run_all.sh" ]; then rm PopLDdecay_run_all.sh ; fi;  for i in "AME" "ASIA" "EUR" "IND"; do for j in $(seq -w 01 20);do k="ChrS$j"; echo "PopLDdecay -MaxDist 100 -MAF 0.05 -Miss 0.2 -i $ld_dir/split_pop/$i/$k.vcf -o $ld_dir/split_pop/$i/$i.$k.stat.gz 2> $ld_dir/split_pop/$i/$i.$k.ld.log" >> PopLDdecay_run_all.sh ; done ; done

sh PopLDdecay_run_all.sh ; rm PopLDdecay_run_all.sh

#to visualize the results from PopLDdecay (codes was modified from https://github.com/BGI-shenzhen/PopLDdecay)
#ml lang/R/4.2.0-foss-2021b #load module @ mogon
for i in "AME" "ASIA" "EUR" "IND"; do  ls $ld_dir/split_pop/$i/*Chr*.stat.gz > $ld_dir/split_pop/$i/$i.chr.list ; Plot_OnePop.pl -inList $ld_dir/split_pop/$i/$i.chr.list -output $ld_dir/split_pop/$i/$i.cat ; done

#per population
for i in "AME" "ASIA" "EUR" "IND"; do echo -e "$ld_dir/split_pop/$i/$i.cat.bin.gz\t$i" ; done > $ld_dir/multi.list

#put all populations all together
Plot_MultiPop.pl -inList $ld_dir/multi.list  -output $ld_dir/4pop.res  -keepR

#put all populations all together, modify the color
cp $ld_dir/4pop.res.r $ld_dir/4pop.res.fix_col.r
sed -i 's/"red"/"#F7931E"/'g $ld_dir/4pop.res.fix_col.r ; sed -i 's/"black"/"#FF0000"/'g $ld_dir/4pop.res.fix_col.r ; sed -i 's/"blue"/"#0000FF"/'g $ld_dir/4pop.res.fix_col.r ; sed -i 's/"Purple"/"#39B54A"/'g $ld_dir/4pop.res.fix_col.r

#regenerate the pdf (png)
Rscript $ld_dir/4pop.res.fix_col.r

mkdir -p ./output; cp $ld_dir/4pop.res.pdf ./output
