#! /usr/bin/perl

use strict;
use warnings;

# my $path=`pwd`;
# chomp $path;

#die "file? perl $0 blast_output6 blast_output0\n" if (@ARGV==0);
die "file? perl $0 blast_output0 original_vcf\n" if (@ARGV==0);

my $min_len=50; #setting minimal aligned sequence length

open OUT0, "<$ARGV[0]";
# Query= ChrS01_13311
# Length=301
# ***** No hits found *****
# Lambda      K        H
#     1.33    0.621     1.12
# Gapped
# Lambda      K        H
#     1.28    0.460    0.850
# Effective search space used: 656737064256
# Query= ChrS01_13356
# Length=301
#                                                                       Score        E
# Sequences producing significant alignments:                          (Bits)     Value
# LG03                                                                  60.2       5e-07
# >LG03
# Length=187626166
#  Score = 60.2 bits (32),  Expect = 5e-07
#  Identities = 69/86 (80%), Gaps = 6/86 (7%)
#  Strand=Plus/Minus
# Query  218       GAACAACGTAAATCCGTGTGTCTAGTTT-TAGAGCGAATTGA--CT-TTTCTTCCCTCTC  273
#                  |||| |||||||||| ||||| || | | | | | | |||||  ||  ||||||||||||
# Sbjct  13993150  GAACCACGTAAATCCTTGTGT-TATTATCTTGTGTG-ATTGATCCTCATTCTTCCCTCTC  13993093
# Query  274       TTTTCTTCTTGTCGTGCAAATTAAGA  299
#                  |||||||||||||||||  || ||||
# Sbjct  13993092  TTTTCTTCTTGTCGTGCGGATCAAGA  13993067
# Lambda      K        H
#     1.33    0.621     1.12
# Gapped
# Lambda      K        H
#     1.28    0.460    0.850
# Effective search space used: 656737064256


my $key;
my $ref_chr;

my $bit;
my $E;
my $ID;
#my %blast_res;
open TMP, ">$ARGV[0].simplified.4ck.tmp";
while (<OUT0>){
        chomp;
        next if (/^\s*$/);

        if (/^Query\=\s*(\w+)$/){
                $key = $1;
                next;
        }
        if (/\*+/){
                $key = undef;
                next;
        }
        if (/^>([\w\.]+)$/){
                $ref_chr=$1;
                #$blast_res{$key}{'chr'}=$ref_chr;
                next;
        }
        if (/^\s*Score\s*\=\s*([^\s]+)\s*bits[^\,]+\,\s*Expect\s*\=\s*([^\s]+)$/){
                $bit=$1;
                $E=$2;
                next;
        }
        if (/^\s*Identities\s+\=\s+([^\s]+)/){
                $ID=$1;
                next;
        }
        if (/^(Query)\s+\d+/){
                #next if (! defined $key);

                #print "$key\n";
                my @a=split/\s+/,$_;

                if ($a[1] <= 151 and $a[3] >= 151 ){
                        #print TMP ">$key\t$ref_chr\t$bit\t$E\t$ID\n";
                        print TMP ">$key\t$ref_chr\t$bit\t$ID\n";
                        print TMP $_,"\n";
                        my $tmp1 = <OUT0>;
                        chomp $tmp1;
                        my $tmp2 = <OUT0>;
                        chomp $tmp2;
                        print TMP "$tmp1\n$tmp2\n";
                }

        }
}
close OUT0;
close TMP;

# Output:
# >ChrS01_13457   LG03    60.2   69/86
# Query  117       GAACAACGTAAATCCGTGTGTCTAGTTT-TAGAGCGAATTGA--CT-TTTCTTCCCTCTC  172
#                  |||| |||||||||| ||||| || | | | | | | |||||  ||  ||||||||||||
# Sbjct  13993150  GAACCACGTAAATCCTTGTGT-TATTATCTTGTGTG-ATTGATCCTCATTCTTCCCTCTC  13993093
# >ChrS01_13487   LG03    60.2   69/ 
# Query  143       TTTTCTTCTTGTCGTGCAAATTAAGA  168
#                  |||||||||||||||||  || ||||
# Sbjct  13993092  TTTTCTTCTTGTCGTGCGGATCAAGA  13993067


open TMP, "<$ARGV[0].simplified.4ck.tmp";
open OUT, ">$ARGV[0].genotype.4ck.tsv";
print OUT "q_chr\tq_pos\ts_chr\ts_start\ts_end\tq_genotype\ts_genotype\tbit_score\tidentity\n";
while (<TMP>){
        chomp;
        if (/^\>(.*)$/){
                my $name = $1;
                $name=~s/_/\t/;
                my @c=split/\t+/,$name;

                my $q = <TMP>;
                my $tmp = <TMP>;
                my $s = <TMP>;
                chomp $q;
                chomp $tmp;
                chomp $s;


                my @a=split/\s+/,$q;
                my @qseq =split//, $a[2];

                my @b=split/\s+/,$s;
                my @sseq =split//, $b[2];

                die "$name not match!\nq:$q\ns:$s\n" if ($q!~/^Query/ || $s!~/^Sbjct/ || @a != 4 || @b != 4);

                my $n = 151- ($a[1] -1);
                my $count_seq=0;
                my $count_gap=0;
                my $loop=$n;

                A: while (0 < 1){
                        for my $i (0..($loop-1)){
                                if ($qseq[$i]=~/[^\-]/){
                                        $count_seq ++;
                                }else {
                                        $count_gap++;
                                }
                        }
                        if ($count_seq < $n){
                                $loop=$count_gap+$n;
                                $count_gap=0;
                                $count_seq=0;
                        } else {
                                last A;
                        }
                }
                my $q_genotype=$qseq[$loop-1];

                my $s_genotype=$sseq[$loop-1];
                #print "$.\tq:$q_genotype\ts:$s_genotype\t$loop\n";
                print OUT "$c[0]\t$c[1]\t$c[2]\t$b[1]\t$b[3]\t$q_genotype\t$s_genotype\t$c[3]\t$c[4]\n";
        }

}
close TMP;
close OUT;

my %hash;

open TSV, "<$ARGV[0].genotype.4ck.tsv" || die "$0; could not open $ARGV[0].genotype.4ck.tsv\n";
# q_chr   q_pos   s_chr   s_start s_end   q_genotype      s_genotype      bit_score       identity
# ChrS01  13457   LG03    13993150        13993093        G       G       60.2    69/86
# ChrS01  13487   LG03    13993092        13993067        T       T       60.2    69/86
# ChrS01  13494   LG03    13993092        13993067        G       G       60.2    69/86
# ChrS01  55994   LG08    172751063       172751122       T       T       222     171/196
# ChrS01  55994   LG13    14448814        14448867        T       T       89.8    52/54
# ChrS01  59763   LG05    18545986        18545927        G       G       206     151/171

while (<TSV>){
        chomp;
        # my @head;
        # if (/^q_chr/){
        #       my @head=split/\t+/,$_;
        #       next;
        # }
        next if (/^q_chr/);
        my @a=split/\t+/,$_;
        my $len;
        if ($a[8]=~/\d+\/(\d+)/){
                $len=$1;
                #print "$len\n";
        } else {
                die "$.: $_ wrong format.\n";
        }
        next if ($a[6] eq "-" || $len < $min_len);
        if (exists $hash{$a[0]}{$a[1]}){
                #print "$hash{$a[0]}{$a[1]}[2]\n";
                if ($a[7] > $hash{$a[0]}{$a[1]}[2]){
                        my @b=($a[5],$a[6],$a[7]);
                        $hash{$a[0]}{$a[1]}=\@b;
                } else {
                        next;
                }
        } else {
                my @b=($a[5],$a[6],$a[7]);
                $hash{$a[0]}{$a[1]}=\@b;
        }
}
close TSV;


open VCF, "<$ARGV[1]";
# ChrS01  18      .       T       A       73070.5 PASS    AC=356;AF=0.813;AN=438;AS_BaseQRankSum=1.600;AS_FS
# ChrS01  53      .       G       T       110720  PASS    AC=379;AF=0.831;AN=456;AS_BaseQRankSum=0.700;AS_FS
open OUTVCF, ">$ARGV[0].4outgroup.vcf";
open OUTTSV, ">$ARGV[0].final.4outgroup.tsv";
print OUTTSV "q_chr\tq_pos\tq_genotype\ts_genotype\tbit_score\n";

open OUTCOUNT, ">$ARGV[0].count.tsv";

open OUTFA, ">$ARGV[0].outgroup.fa";
print OUTFA ">outgroup\n";

my $vcf_sites=0; #num of SNP sites in original VCF
my $found_sites=0; #num of orthologous sites found in outgroup sp.
my $identical_sites=0; #num of orthologous sites found in outgroup sp. that indentical to targeted sp.

while (<VCF>){
        chomp;
        if (/^#/){
                print OUTVCF "$_\n";
                next;
        }
        $vcf_sites ++;
        my @a=split/\t/,$_;
        if (exists $hash{$a[0]}{$a[1]}){
                $found_sites++;
                my $char=$hash{$a[0]}{$a[1]}[0];
                $char=uc($char);
                $a[3]=uc($a[3]);
                die "$_ not match.\n$a[0]\t$a[1]\t$char\n" if ($char ne $a[3]);
                print OUTVCF "$_\n";
                print OUTFA "$hash{$a[0]}{$a[1]}[1]";
                
                my @b=@{$hash{$a[0]}{$a[1]}};
                $b[0]=uc($b[0]);
                $b[1]=uc($b[1]);
                
                if ($b[0] eq $b[1]){
                        $identical_sites++;
                }
                my $tmp =join("\t",@b);
                print OUTTSV "$a[0]\t$a[1]\t$tmp\n";
        }
}
print  OUTFA "\n";

print OUTCOUNT "Num of mutations in VCF:\t$vcf_sites\n";
print OUTCOUNT "Num of mutations found in outgroup sp.:\t$found_sites\n";
print OUTCOUNT "Num of identical sites:\t$identical_sites\n";
my $p = $found_sites / $vcf_sites * 100;
print OUTCOUNT "Perc. of mutations found in outgroup sp.:\t$p\n";
$p = $identical_sites / $found_sites * 100;
print OUTCOUNT "Perc. of identical sites in all found:\t$p\n";
$p = $identical_sites / $vcf_sites * 100;
print OUTCOUNT "Perc. of identical sites in all sites:\t$p\n";

close OUTCOUNT;
close VCF;
close OUTVCF;
close OUTFA;
close OUTTSV;


__END__
explain of blast output6 format:
 1.      qseqid  query (e.g., unknown gene) sequence id
 2.      sseqid  subject (e.g., reference genome) sequence id
 3.      pident  percentage of identical matches
 4.      length  alignment length (sequence overlap)
 5.      mismatch        number of mismatches
 6.      gapopen         number of gap openings
 7.      qstart  start of alignment in query
 8.      qend    end of alignment in query
 9.      sstart  start of alignment in subject
 10.     send    end of alignment in subject
 11.     evalue  expect value
 12.     bitscore        bit score
