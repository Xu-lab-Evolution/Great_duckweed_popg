#! /usr/bin/perl

use strict;
use warnings;



#die "file? perl $0 blast_output6 blast_output0\n" if (@ARGV==0);
die "snp_list, vcf & /path/to/vcf2phylip.py?\n" if (@ARGV!=3);

die "$ARGV[0] name not match! please provide the snp_list in relative path (./*final.4outgroup.tsv.snp.list) or abs path\n" unless ($ARGV[0]=~/^.*\/(\w+)\.[^\/]+$/);
my $run_name =$1;
chomp $run_name;
$run_name=$run_name.".run_raxml";

my $path=`dirname $ARGV[0]`;
chomp $path;
print "$path\n";

my $vcf2phylip_path=$ARGV[2];
chomp $vcf2phylip_path;

my $output_dir="$path/$run_name";
`mkdir -p  $output_dir`;
#chdir "$path/$run_name";

my %hash;
open IN, "<$ARGV[0]";
while (<IN>){
        chomp;
        $hash{$_}=0;
}
close IN;

my @ck;
open VCF, "<$ARGV[1]";
open OUTVCF, ">$output_dir/$run_name.tmp.vcf";
while (<VCF>){
        chomp;
        if (/^#/){
                print OUTVCF "$_\n";
                next;
        }
        my @a=split/\t/,$_;
        if (exists $hash{$a[0]."_".$a[1]}){
                print OUTVCF "$_\n";
                push @ck, $a[0]."_".$a[1];
        }
}
close VCF;
close OUTVCF;

print "$vcf2phylip_path/vcf2phylip.py  -f -n  -i  $output_dir/$run_name.tmp.vcf --output-folder $output_dir\n";
`$vcf2phylip_path/vcf2phylip.py  -f -n  -i  $output_dir/$run_name.tmp.vcf --output-folder $output_dir`;

open TSV_LIST, "<$ARGV[0].file_list";
my $c=0;
while (<TSV_LIST>){
        chomp;
        open IN, "<$_";
        $c++;
        my $sp_name;
        die "$_ name not match!\n" unless (/^.*\/(\w+)\.[^\/]+$/);
        $sp_name=$1;
        open FA, ">$output_dir/$c.$sp_name.tmp.fa";
        print FA ">$sp_name\n";
# q_chr   q_pos   q_genotype      s_genotype      bit_score
# ChrS01  51356   C       T       252
# ChrS01  51389   C       C       202
        while (<IN>){
                chomp;
                next if (/^q_chr/);
                my @a=split/\t/,$_;
                if (exists $hash{$a[0]."_".$a[1]}){
                        print FA "$a[3]";
                }
        }
        print FA "\n";
        close FA;
        close IN;
}
close TSV_LIST;

`cat $output_dir/*.tmp.fa $output_dir/$run_name.tmp.min4.fasta > $output_dir/$run_name.fa`;

open MDTEST, ">$run_name.MDtest.sh";
print MDTEST "#! /bin/bash

modeltest-ng -i -p 8 $output_dir/$run_name.fa -o $output_dir/$run_name.MDtest\n";

close MDTEST;

#GTR+I+G4 model was selected based on the modeltest's results
open RAXML, ">$run_name.RAxML.sh"; 
print RAXML "#! /bin/bash

raxml-ng --msa $output_dir/$run_name.fa --model TVM+G4  --all  --threads 40 --bs-metric fbp,tbe  --seed 142857 --tree pars{50},rand{50} --prefix $output_dir/$run_name.raxml\n";

close RAXML;

__END__

