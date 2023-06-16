#! /bin/bash

#ml bio/BLAST+/2.12.0-gompi-2021b  #@ mogon cluster, change it when naccesary

#cd ../../../data//Fig_1/01.Phylogenetic_tree//
makeblastdb  -in ../../../data//Fig_1/01.Phylogenetic_tree//Colocasia_esculenta.Taro_Lachesis_assembly_Chr.fa -dbtype nucl  -max_file_sz 3GB  -logfile ../../../data//Fig_1/01.Phylogenetic_tree//Colocasia_esculenta.Taro_Lachesis_assembly_Chr.fa.log

#cd ../../../data//Fig_1/01.Phylogenetic_tree//
blastn -db ../../../data//Fig_1/01.Phylogenetic_tree//Colocasia_esculenta.Taro_Lachesis_assembly_Chr.fa  -max_target_seqs 1 -query ../../../data//Fig_1/01.Phylogenetic_tree///SP_228.basic_set.snp.recode.rm_cluster3-10.vcf.PASS.vcf.flanking150.fa -outfmt 0  -evalue 1e-6  -num_threads 40  -out ../../../data//Fig_1/01.Phylogenetic_tree///Colocasia_esculenta.flanking150.blast.fmt0.res

