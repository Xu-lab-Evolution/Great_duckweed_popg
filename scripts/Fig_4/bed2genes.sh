bedtools=/home/pablo/Programs/bedtools2/bin/bedtools
ann=/home/pablo/Documents/Duckweed/annotation/SpGA2022.3_nopseudo0_onlyGenes.gff3
#dir=outputs_EUR
dir=outputs_IND

top01percent_anc=${dir}/results_ancestral.bed
#top01percent_eur=${dir}/results_europe.bed
top01percent_ind=${dir}/results_india.bed
top01percent_asi=${dir}/results_asia.bed

$bedtools intersect -a $ann -b $top01percent_anc > ${dir}/out_3PCLR_genes_ancestral.gff
#$bedtools intersect -a $ann -b $top01percent_eur > ${dir}/out_3PCLR_genes_europe.gff
$bedtools intersect -a $ann -b $top01percent_ind > ${dir}/out_3PCLR_genes_india.gff
$bedtools intersect -a $ann -b $top01percent_asi > ${dir}/out_3PCLR_genes_asia.gff


#pop=("ancestral" "europe" "asia")
pop=("ancestral" "india" "asia")
echo -e "Gene ID\tChromosome\tStart\tEnd\tFunction" > header.txt

for str in ${pop[@]}
do
	infile=${dir}/out_3PCLR_genes_$str.gff
	cut -f 9 $infile | cut -f 1 -d ";" | sed 's/ID=//g' > column1.txt
	cut -f 1,4,5 $infile > column2-4.txt
	cut -f 9 $infile | cut -f 3 -d ";" | awk '{ if ($0 ~ /Note=/) {print $0} else {print "-"} }' | sed 's/Note=//g' > column5.txt
	paste column1.txt column2-4.txt column5.txt > columns.txt
	cat header.txt columns.txt > ${dir}/table_genes1percent_3PCLR_$str.txt
	rm column*
done

