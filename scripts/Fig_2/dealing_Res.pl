#! /usr/bin/perl

use strict;
use warnings;
use List::MoreUtils qw(zip);

die "file?\n" if (@ARGV==0);

my %hash;
my %genome;
my %genome_sum;
my $dir=`dirname $ARGV[0]`;
chomp $dir;


open LIST, "<$ARGV[0]";
open ERR, ">$dir/completeness.ck";
my @name;

my %cor=("sites","pi","sites_coding","pi_coding","sites_noncoding","pi_noncoding","N_sites","piN","S_sites","piS");

while (<LIST>){
	chomp;
	open IN, "<$_";
	while (<IN>){
		chomp;
		my @a=split/\t/,$_;
		my $tmp=@a+1;
		die "file format wrong: $_\n$tmp elements exist, not 22\n" if (@a != 22);
		
		my @value;
		my %tmp;
		my $pop;
		my $chr;
		if ($a[0]=~m/^file/){
			@name=@a;
			next;
		} else {
			#die "file cotent wrong: $_\n" if ($a[4]!~m/^[\d\.]+$/);
			my @b = split/\//, $a[0];
			$pop=$b[-3];
			$chr=$b[-2];

			#print "ck: $pop\t$chr\n";
			@value=@a;
			print ERR "$pop/$chr\t$a[4]\n";
		}
		#%tmp=zip @name,@value;
		#@tmp{ @name } = @value ;
		for ( my $i = 0 ; $i <= $#name; $i++) {
			my $k = $name[$i];
			my $v = $value[$i];
			$tmp{$k} =$v;
		}
		#print "ck: $pop\t$chr\n";
		$hash{$pop}{$chr}{"pi"} = $tmp{"pi"};
		$hash{$pop}{$chr}{"pi_coding"} = $tmp{"pi_coding"};
		$hash{$pop}{$chr}{"pi_noncoding"} = $tmp{"pi_noncoding"};
		$hash{$pop}{$chr}{"piN"} = $tmp{"piN"};
		$hash{$pop}{$chr}{"piS"} = $tmp{"piS"};

		$genome{$chr}{$cor{"sites_coding"}} = $tmp{"sites_coding"};
		$genome{$chr}{$cor{"sites_noncoding"}} = $tmp{"sites_noncoding"};
		$genome{$chr}{$cor{"sites"}} = $tmp{"sites"};
		$genome{$chr}{$cor{"N_sites"}} = $tmp{"N_sites"};
		$genome{$chr}{$cor{"S_sites"}}= $tmp{"S_sites"};
	}
}

close LIST;
close ERR;

my @keys=("pi","pi_coding","pi_noncoding","piN","piS");
my @pop=("AME","ASIA","EUR","IND","ALLPOP");
my @chr=("ChrS01","ChrS02","ChrS03","ChrS04","ChrS05","ChrS06","ChrS07","ChrS08","ChrS09","ChrS10","ChrS11","ChrS12","ChrS13","ChrS14","ChrS15","ChrS16","ChrS17","ChrS18","ChrS19","ChrS20");

for my $c (@chr){
	for my $k (keys %{$genome{$c}}){
		if (not exists $genome_sum{$k}){
			$genome_sum{$k}=0;
		}

		$genome_sum{$k}+=$genome{$c}{$k};
	}
}

open GENOME_STAT, ">$dir/Genome_check.stat";
for my $k (keys %genome_sum){
	print GENOME_STAT "$k\t$genome_sum{$k}\n";
}
print GENOME_STAT "\n";
my @kk=keys %genome_sum;
print GENOME_STAT (join"\t",@kk)."\n";

for my $c (@chr){
	print GENOME_STAT "$c";
	for my $k (@kk){
		print GENOME_STAT "\t$genome{$c}{$k}";
	}
	print GENOME_STAT "\n";
}
close GENOME_STAT;



my %sum;


my $id = join "\t",@keys;

open ERR, ">$dir/err.ck";
for my $p (@pop){
	open TMP, ">$dir/SNPGenie_res_dealing.tmp.$p";
	
	print TMP "Chr\t$id\n";

	A: for my $c (@chr){
		if (! exists $hash{$p}{$c}){
			print ERR "$p/$c\n";
			next A;
		} 
		print TMP "$c";
		for my $k (@keys){
			#print "ck: $p\t$c\t$k\n";

			print TMP "\t$hash{$p}{$c}{$k}";
			$sum{$p}{$k}+=$hash{$p}{$c}{$k} * $genome{$c}{$k};
		}
		print TMP "\n";
	}
	close TMP;
}
close ERR;

open OUT, ">$dir/SNPGenie_res_dealing.res";
print OUT "Pop\t$id\n";
for my $p (@pop){
	print OUT "$p";
	for my $k (@keys){
		my $r;
		if ($genome_sum{$k}==0){
			$r=0;
		} else {
		 
			$r=$sum{$p}{$k} / $genome_sum{$k};
		}
		print OUT "\t$r";
	}
	print OUT "\n";
}
close OUT;


__END__
file	sites	sites_coding	sites_noncoding	pi	pi_coding	pi_noncoding	N_sites	S_sites	piN	piS	mean_dN_vs_ref	mean_dS_vs_ref	mean_gdiv_polymorphic	mean_N_gdiv	mean_S_gdiv	mean_gdiv	sites_polymorphic	mean_gdiv_coding_poly	sites_coding_poly	mean_gdiv_noncoding_poly	sites_noncoding_poly
/scratch/tmp/ywang1/01.duckweed/popgenomics/neutral_stat/split_chr//AME/ChrS19/ChrS19.vcf.AF_corrected.vcf	3727809	473454	3254355	0.000974306450639792	0.0006604793088573	0.00101996302946586	152140.280237587	49129.334599811	0.000353927851349201	0.000858030951549198	0.000309208666744975	0.000714132564777546	0.00660900161906613	0.26518144702167	0.26512080390582	0.000961614664760828	14345	0.259463153635798	1190	0.249027339338719	13155
