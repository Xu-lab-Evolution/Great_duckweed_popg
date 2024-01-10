#! /usr/bin/perl

use strict;
use warnings;

die "Usage: perl $0 vcf ref_dir gtf_dir snpgenie_dir switch\n" if (@ARGV < 4);

my $vcf_file=$ARGV[0];
my $vcf_dir=`dirname $vcf_file`;
my $chr=`basename  $vcf_file .vcf`;

my $ref_dir=$ARGV[1];
my $gtf_dir=$ARGV[2];
my $snpg_dir=$ARGV[3];
my $switch=$ARGV[4];

chomp($vcf_file,$vcf_dir,$ref_dir,$gtf_dir,$snpg_dir,$chr,$switch);

#print "ref: $ref_dir\ngtf: $gtf_dir\n";

if ($switch=~m/\w/){ #swith of whether to generate sbatch cmd line on mogon
        open CMD, ">$vcf_dir/sbatch.snpgenie.sh";
        print CMD "\#\!/bin/bash
#SNPgenie

\#SBATCH --export=ALL              
\#SBATCH -A m2_jgu-EvolTroph
\#SBATCH --partition parallel
\#SBATCH --nodes=1                  
\#SBATCH -c 1        
\#SBATCH --mem=12G  
\#SBATCH --time=01:00:00             # the max wallclock time (time limit your job will run)
\#SBATCH --job-name=SNPgenie.$chr
\#SBATCH --output=$vcf_dir/SNPgenie.out
\#SBATCH --error=$vcf_dir/SNPgenie.err

ml lang/Perl/5.36.0-GCCcore-12.2.0

perl $snpg_dir/correct_AF_info.4vcf.pl $vcf_file
perl $snpg_dir/snpgenie.pl  --vcfformat=1 --snpreport=$vcf_file.AF_corrected.vcf --fastafile=$ref_dir/$chr.fa --gtffile=$gtf_dir/$chr.gtf

ml bio/VCFtools/0.1.16-GCC-11.2.0
vcftools --vcf $vcf_file.AF_corrected.vcf --site-pi --out $vcf_file.vcftools.pi.stat\n";
        close CMD;

        open ALL, ">>./snpgenie.all_run.list.sh";
        print ALL "cd $vcf_dir; sbatch $vcf_dir/sbatch.snpgenie.sh\n";
        close ALL;  
} else {
        open CMD, ">$vcf_dir/snpgenie.sh";
        print CMD "\#\!/bin/bash

#ml lang/Perl/5.36.0-GCCcore-12.2.0

perl correct_AF_info.4vcf.pl $vcf_file
perl $snpg_dir/snpgenie.pl  --vcfformat=1 --snpreport=$vcf_file.AF_corrected.vcf --fastafile=$ref_dir/$chr.fa --gtffile=$gtf_dir/$chr.gtf

#ml bio/VCFtools/0.1.16-GCC-11.2.0
vcftools --vcf $vcf_file --site-pi --out $vcf_file.vcftools.pi.stat\n";
        close CMD;

        open ALL, ">>./snpgenie.all_run.list.sh";
        print ALL "cd $vcf_dir; sh $vcf_dir/snpgenie.sh\n";
        close ALL;
}




