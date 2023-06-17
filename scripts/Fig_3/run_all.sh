data_dir="../../data/" #diractory stores input VCF and temporary VCF
general_scripts_dir="../other_scripts/"
DNA_methylation_dir="$data_dir/Fig_3/00.calculate_wML/"


######################################perform binormial test filtration
#lang/R/4.2.0-foss-2021b #@ mogon
for i in `cat $DNA_methylation_dir/samp.list`; do  zcat $DNA_methylation_dir/CX_report/${i}_CX_report.txt.gz | grep -v "JN160603.2\|JQ804980.1\|pseudo0" | gzip > $DNA_methylation_dir/CX_report/${i}_CX_report.txt.nuc.gz ; done

for i in `cat $DNA_methylation_dir/samp.list`; do Rscript Binomial_Test.filtration.4Bismark_cov.R $DNA_methylation_dir/convs.20samp.tab $DNA_methylation_dir/CX_report/${i}_CX_report.txt.nuc.gz ; done

for i in `cat $DNA_methylation_dir/samp.list`; do gunzip -c $DNA_methylation_dir/CX_report/${i}_CX_report.filtered.gz >  $DNA_methylation_dir/CX_report/${i}_CX_report.filtered; done


######################################calculate weighted methylation level for 1)genome-wide 2)gene region 3)TE region; calculate the proportion of methylated cytosine

mkdir -p $DNA_methylation_dir/calc_ML_res; for i in `cat $DNA_methylation_dir/samp.list`; do Rscript  calc_regional_wML.prop_mC.R  $DNA_methylation_dir/CX_report/${i}_CX_report.filtered $data/Spolyrhiza_annotation_SPGA2022/4_methyl $DNA_methylation_dir/calc_ML_res ; done

######################################for global weighted methylation level visualization
Rscript global_ML.barplot.R $DNA_methylation_dir/samp.list $DNA_methylation_dir/calc_ML_res/regional_wML_tsv

######################################to generate meta plots of methylation level over genes and TEs

obj="Gene"
for cg_context in "CG" "CHG" "CHH" ; do ViewBS MethOverRegion --region $data/Spolyrhiza_annotation_SPGA2022/4_methyl/SpGA2022.gene.bed --sample $DNA_methylation_dir//CX_report/AM11_CX_report.filtered.gz,AM11   --sample $DNA_methylation_dir//CX_report/AM32_CX_report.filtered.gz,AM32   --sample $DNA_methylation_dir//CX_report/AM4_CX_report.filtered.gz,AM4   --sample $DNA_methylation_dir//CX_report/AM85_CX_report.filtered.gz,AM85   --sample $DNA_methylation_dir//CX_report/AM87_CX_report.filtered.gz,AM87 --sample $DNA_methylation_dir//CX_report/AS115_CX_report.filtered.gz,AS115   --sample $DNA_methylation_dir//CX_report/AS150_CX_report.filtered.gz,AS150   --sample $DNA_methylation_dir//CX_report/AS165_CX_report.filtered.gz,AS165   --sample $DNA_methylation_dir//CX_report/AS182_CX_report.filtered.gz,AS182   --sample $DNA_methylation_dir//CX_report/AS72_CX_report.filtered.gz,AS72   --sample $DNA_methylation_dir//CX_report/EU15_CX_report.filtered.gz,EU15   --sample $DNA_methylation_dir//CX_report/EU23_CX_report.filtered.gz,EU23   --sample $DNA_methylation_dir//CX_report/EU31_CX_report.filtered.gz,EU31   --sample $DNA_methylation_dir//CX_report/EU47_CX_report.filtered.gz,EU47   --sample $DNA_methylation_dir//CX_report/EU55_CX_report.filtered.gz,EU55   --sample $DNA_methylation_dir//CX_report/IN227_CX_report.filtered.gz,IN227 --sample $DNA_methylation_dir//CX_report/IN236_CX_report.filtered.gz,IN236   --sample $DNA_methylation_dir//CX_report/IN39_CX_report.filtered.gz,IN39   --sample $DNA_methylation_dir//CX_report/IN43_CX_report.filtered.gz,IN43   --sample $DNA_methylation_dir//CX_report/IN59_CX_report.filtered.gz,IN59  --outdir $DNA_methylation_dir/ViewBS/methOver.${obj} --prefix 20samp.${obj}.$cg_context --context $cg_context --regionName $obj 2> $DNA_methylation_dir/ViewBS/4.MethOverRegion.${obj}.${cg_context}.out ; done

obj="TE"
for cg_context in "CG" "CHG" "CHH" ; do ViewBS MethOverRegion --region $data/Spolyrhiza_annotation_SPGA2022/4_methyl/RE.bed4.f --sample $DNA_methylation_dir//CX_report/AM11_CX_report.filtered.gz,AM11   --sample $DNA_methylation_dir//CX_report/AM32_CX_report.filtered.gz,AM32   --sample $DNA_methylation_dir//CX_report/AM4_CX_report.filtered.gz,AM4   --sample $DNA_methylation_dir//CX_report/AM85_CX_report.filtered.gz,AM85   --sample $DNA_methylation_dir//CX_report/AM87_CX_report.filtered.gz,AM87 --sample $DNA_methylation_dir//CX_report/AS115_CX_report.filtered.gz,AS115   --sample $DNA_methylation_dir//CX_report/AS150_CX_report.filtered.gz,AS150   --sample $DNA_methylation_dir//CX_report/AS165_CX_report.filtered.gz,AS165   --sample $DNA_methylation_dir//CX_report/AS182_CX_report.filtered.gz,AS182   --sample $DNA_methylation_dir//CX_report/AS72_CX_report.filtered.gz,AS72   --sample $DNA_methylation_dir//CX_report/EU15_CX_report.filtered.gz,EU15   --sample $DNA_methylation_dir//CX_report/EU23_CX_report.filtered.gz,EU23   --sample $DNA_methylation_dir//CX_report/EU31_CX_report.filtered.gz,EU31   --sample $DNA_methylation_dir//CX_report/EU47_CX_report.filtered.gz,EU47   --sample $DNA_methylation_dir//CX_report/EU55_CX_report.filtered.gz,EU55   --sample $DNA_methylation_dir//CX_report/IN227_CX_report.filtered.gz,IN227 --sample $DNA_methylation_dir//CX_report/IN236_CX_report.filtered.gz,IN236   --sample $DNA_methylation_dir//CX_report/IN39_CX_report.filtered.gz,IN39   --sample $DNA_methylation_dir//CX_report/IN43_CX_report.filtered.gz,IN43   --sample $DNA_methylation_dir//CX_report/IN59_CX_report.filtered.gz,IN59  --outdir $DNA_methylation_dir/ViewBS/methOver.${obj} --prefix 20samp.${obj}.$cg_context --context $cg_context --regionName $obj 2> $DNA_methylation_dir/ViewBS/4.MethOverRegion.${obj}.${cg_context}.out ; done

####to use ggplot2 to draw the meta plots
Rscript 20samp.meta_plot.TE_Gene.R $DNA_methylation_dir/ViewBS/ ./output
