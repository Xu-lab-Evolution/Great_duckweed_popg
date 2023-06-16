#! /usr/bin/perl

use strict;
use warnings;

#my $path=`pwd`;
#chomp $path;

die "file? perl $0 /path/to/out_group_ref.fa /path/to/original_snp.vcf /path/to/ref.fa /path/to/output \n" if (@ARGV!=4);
#die "$ARGV[1] name or path not match (please specify absolute path)!\n" unless $ARGV[1]=~/^(\/.*\/)\w+\.[^\/]+$/;
 
#my $vcf_path=`dirname $ARGV[1]`;

#die "$ARGV[0] name or path not match (please specify absolute path)!\n" unless $ARGV[0]=~/^(\/.*\/)(\w+)\.[^\/]+$/;
#$ARGV[0]=~/^(\/.*\/)(\w+)\.[^\/]+$/;
#my $out_group_ref_path=`dirname $ARGV[0]`;
my $outgroup_name=`basename $ARGV[0]`;

$outgroup_name=~s/\..*$//;
chomp $outgroup_name;

my $out_group_ref_path=$ARGV[3];
#print "dealing with outgroup (name: $outgroup_name)\n";

print "Outgroup species: $outgroup_name\n";

# my $path="$out_dir/$outgroup_name";
# `mkdir -p $path`;

my $vcf_file=`basename $ARGV[1]`;
chomp $vcf_file;

my $flk_seq="$out_group_ref_path/$vcf_file.flanking150.fa";
if (! -e $flk_seq){
        print "$flk_seq not found, try to creat...\n";

        #my $i=$vcf_file;

        open TMP, ">./get_flanking150.sh";
        print TMP "#module load GCC/8.2.0-2.31.1 BEDTools/2.28.0
#ml bio/SAMtools/1.14-GCC-11.2.0 #\@ mogon cluster, change it when naccesary
samtools faidx $ARGV[2]
awk '{OFS=\"\\t\"; print \$1,\$2}' $ARGV[2].fai > $ARGV[2].genome.txt
awk '{OFS=\"\\t\"; if (!/^#/){print \$1,\$2-1,\$2,\$1.\"_\"\$2,\$4\"/\"\$5}}' $ARGV[1] > $out_group_ref_path/$vcf_file.bed

#ml bio/BEDTools/2.29.2-iccifort-2020.1.217 #\@ mogon cluster, change it when naccesary
bedtools slop -i $out_group_ref_path/$vcf_file.bed  -g $ARGV[2].genome.txt -l 150 -r 150  > $out_group_ref_path/$vcf_file.flanking150.bed
awk '{if(\$3-\$2==301){print}}' $out_group_ref_path/$vcf_file.flanking150.bed | cat > $out_group_ref_path/$vcf_file.flanking150.fixed.bed
bedtools getfasta  -fi $ARGV[2]  -name -bed $out_group_ref_path/$vcf_file.flanking150.fixed.bed > $out_group_ref_path/$vcf_file.flanking150.fa\n";

`sleep 5`;
`sh ./get_flanking150.sh`;
        close TMP;
}



open OUT ,">./$outgroup_name.find_OG.sh";
# print OUT "#! /bin/bash
# #SBATCH --export=ALL               # Start with a clean environment
# #SBATCH --nodes=1                   # the number of nodes you want to reserve
# #SBATCH -c=40        # the number of CPU cores per node
# #SBATCH --mem=20G                 # how much memory is needed per node (units can be: K, M, G, T)
# #SBATCH --partition=normal          # on which partition to submit the job
# #SBATCH --time=4:00:00             # the max wallclock time (time limit your job will run)
# #SBATCH --job-name=$outgroup_name.find_OG        # the name of your job
# #SBATCH --output=$outgroup_name.find_OG.output.dat         # the file where output is written to (stdout & stderr)
# module load GCC/8.2.0-2.31.1  OpenMPI/3.1.3 icc/2019.1.144-GCC-8.2.0-2.31.1  impi/2018.4.274 ifort/2019.1.144-GCC-8.2.0-2.31.1  impi/2018.4.274 BLAST+/2.9.0

# cd $out_group_ref_path
# makeblastdb  -in $ARGV[0] -dbtype nucl  -max_file_sz 3GB   -logfile $ARGV[0].log

# cd $path
# blastn -db $ARGV[0]  -max_target_seqs 1 -query $flk_seq -outfmt 0  -evalue 1e-6  -num_threads 70  -out $path/$outgroup_name.flanking150.blast.fmt0.res
# \n";


print OUT "#! /bin/bash

#ml bio/BLAST+/2.12.0-gompi-2021b  #\@ mogon cluster, change it when naccesary

#cd $out_group_ref_path
makeblastdb  -in $ARGV[0] -dbtype nucl  -max_file_sz 3GB  -logfile $ARGV[0].log

#cd $out_group_ref_path
blastn -db $ARGV[0]  -max_target_seqs 1 -query $flk_seq -outfmt 0  -evalue 1e-6  -num_threads 40  -out $out_group_ref_path/$outgroup_name.flanking150.blast.fmt0.res
\n";

close OUT;

__END__

