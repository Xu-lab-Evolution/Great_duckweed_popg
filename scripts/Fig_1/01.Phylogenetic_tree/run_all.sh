data_dir="../../../data/" #diractory stores input VCF and temporary VCF
scripts_dir=`pwd` #"scripts/Fig_1/01.Phylogenetic_tree/"
general_scripts_dir="../../other_scripts/"
basic_set_SNP_vcf_dcprs="$data_dir/SNP_vcf/SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf"
phyl_dir="$data_dir/Fig_1/01.Phylogenetic_tree"


#this step worked for two purposes: 1) extract 150bp flanking regions in fasta file of each SNP from VCF; 2) prepare blast script, put outgroup referece as target database while SNPs and corresponding flanking stretchs as query sequences. 

perl find_OG_pipeline.integrated.pl $phyl_dir/Colocasia_esculenta.Taro_Lachesis_assembly_Chr.fa $basic_set_SNP_vcf_dcprs  $data_dir/Spolyrhiza_ref/SP_combined.fasta $phyl_dir/

#to blast, The results from blastn should be stored in "outgroup_name.flanking150.blast.fmt0.res"
sh Colocasia_esculenta.find_OG.sh

#to check, and get outgroup fasta file (and also corresponding SNP list)
perl get_othologous_genotype.pl $phyl_dir/Colocasia_esculenta.flanking150.blast.fmt0.res $basic_set_SNP_vcf_dcprs

grep -v "q_chr" $phyl_dir/Colocasia_esculenta.flanking150.blast.fmt0.res.final.4outgroup.tsv | awk '{print $1"_"$2}' > $phyl_dir/Colocasia_esculenta.flanking150.blast.fmt0.res.final.4outgroup.tsv.snp.list

#to extract SNP from VCF file as a fasta, and combine original SNPs fasta with outgroup fasta to form a new one, which is prepared as the input for raxml-ng
echo "$phyl_dir/Colocasia_esculenta.flanking150.blast.fmt0.res.final.4outgroup.tsv" > $phyl_dir/Colocasia_esculenta.flanking150.blast.fmt0.res.final.4outgroup.tsv.snp.list.file_list

perl produce_fa.from_VCF_with_OG.pl $phyl_dir/Colocasia_esculenta.flanking150.blast.fmt0.res.final.4outgroup.tsv.snp.list $basic_set_SNP_vcf_dcprs $general_scripts_dir

#the previous step generated the fasta file which contains all SNPs that could be found from both target species and outgroup species, as well bash scripts to run modeltest-ng and raxml-ng

#run modeltest and raxml

sh Colocasia_esculenta.run_raxml.MDtest.sh

sh Colocasia_esculenta.run_raxml.RAxML.sh

#collect result
mkdir -p ./output; cp $phyl_dir/Colocasia_esculenta.run_raxml/Colocasia_esculenta.raxml.raxml.supportTBE ./output

#the final phylogeny was plotted using ITOL
