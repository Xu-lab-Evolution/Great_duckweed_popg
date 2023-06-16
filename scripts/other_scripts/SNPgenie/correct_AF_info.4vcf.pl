#! /usr/bin/perl

use strict;
use warnings;

die "file?\n" if (@ARGV==0);
open VCF, "<$ARGV[0]";
open OUT, ">$ARGV[0].AF_corrected.vcf";
open LOG, ">$ARGV[0].AF_corrected.log";

A: while (<VCF>){
	chomp;
	if (/^#/){
		print OUT "$_\n";
		next A;
	} else {
		my @a=split/\t/,$_;

		my @info=split/;/,$a[7];
		my $AN="NA";
		my $AF="NA";
		for my $entry (@info){
			if ($entry=~m/AN\=([\d\.]+)/){
				$AN=$1;
			}
			if ($entry=~m/AF\=([\d\.]+)/){
				$AF=$1;
			}
		}

		if ($AN eq "NA" and $AF eq "NA") {
			die "Line: $. dosen't have AF or AN record:\n$_\n";
		}

		my @genotype=@a[9..(@a-1)];
		my $n_allele=0;
		my $n_alt_allele=0;
		B: for my $i (@genotype){
			my @b=split/:/,$i;
			if ($b[0]=~/([\d\.])[\/\|]([\d\.])/ ){
				my $first_allele=$1;
				my $sec_allele=$2;
				if ($first_allele=~/\d/){
					$n_allele += 2;
					if ($first_allele > 0){
						$n_alt_allele += 1;
					}
					if ($sec_allele > 0){
						$n_alt_allele += 1;
					}
				}
			} elsif ($b[0]=~/^[\.\d]$/){
				print LOG "warning: Line: $., $i dosen't match '0(1)/(|)0(1)' regulation\n";
				#next B;
			} else {
				die "Line: $., $i dosen't match '0(1)/(|)0(1)' regulation\n";
			}
		}

		if ($n_allele > 0 and $n_alt_allele > 0){
			#$AF=sprintf("%.8f",$n_alt_allele/$n_allele);
			$AF=$n_alt_allele/$n_allele;
		} else {
			$AF=0;
		}
		
		#$a[7]="AN=$n_allele;AF=$AF";
		my $info=$a[7];
		if ($info=~s/AN=\d+\;//){
			if($info=~s/AF=[\d\-e\.]+\;//) {
				$a[7]="AN=$n_allele;AF=$AF;"."$info";
			} else {
				die "Line: $., AF INFO wrong\n";
			}
		} else {
			die "Line: $., AN INFO wrong\n";
		}

		#$a[7]=$info;
	
		my $out=join "\t", @a;
		print OUT "$out\n";
	}
}

print LOG "\nCorrection completed.\n";

close VCF;
close OUT;
close LOG;
