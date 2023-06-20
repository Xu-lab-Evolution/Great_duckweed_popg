#!/bin/bash

#SBATCH -J sweed                  # Job name
#SBATCH -o sweed.%j.out           # Specify stdout output file (%j expands to jobId)
#SBATCH -p smp                   # Queue name 'smp' or 'parallel' on Mogon II
#SBATCH -n 1                     # Total number of tasks, here explicitly 1
#SBATCH --mem 2G                 # The default is 300M memory per job. You'll likely have to adapt this to your needs
#SBATCH -t 3:00:00              # Run time (hh:mm:ss)
 
#SBATCH -A  m2_jgu-evoltroph     # Specify allocation to charge against

#SBATCH --mail-type=END
#SBATCH --mail-user=pduchnbo@uni-mainz.de
 
# Load all necessary modules if needed (these are examples)
# Loading modules in the script ensures a consistent environment.
 
# Launch the executable
sweed=/lustre/project/m2_jgu-evoltroph/pduchnbo/Programs/sweed-master/SweeD

$sweed -name ChrS01 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS01.159.samp.0missing.vcf -grid 11467
$sweed -name ChrS02 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS02.159.samp.0missing.vcf -grid 8941
$sweed -name ChrS03 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS03.159.samp.0missing.vcf -grid 8796
$sweed -name ChrS04 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS04.159.samp.0missing.vcf -grid 8492
$sweed -name ChrS05 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS05.159.samp.0missing.vcf -grid 8390
$sweed -name ChrS06 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS06.159.samp.0missing.vcf -grid 8131
$sweed -name ChrS07 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS07.159.samp.0missing.vcf -grid 8108
$sweed -name ChrS08 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS08.159.samp.0missing.vcf -grid 7340
$sweed -name ChrS09 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS09.159.samp.0missing.vcf -grid 7208
$sweed -name ChrS10 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS10.159.samp.0missing.vcf -grid 7041
$sweed -name ChrS11 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS11.159.samp.0missing.vcf -grid 6553
$sweed -name ChrS12 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS12.159.samp.0missing.vcf -grid 5946
$sweed -name ChrS13 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS13.159.samp.0missing.vcf -grid 5477
$sweed -name ChrS14 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS14.159.samp.0missing.vcf -grid 5104
$sweed -name ChrS15 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS15.159.samp.0missing.vcf -grid 4726
$sweed -name ChrS16 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS16.159.samp.0missing.vcf -grid 4624
$sweed -name ChrS17 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS17.159.samp.0missing.vcf -grid 4565
$sweed -name ChrS18 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS18.159.samp.0missing.vcf -grid 4370
$sweed -name ChrS19 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS19.159.samp.0missing.vcf -grid 3728
$sweed -name ChrS20 -input /lustre/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/selection_scan_LASSI/CHR_VCF_without_MISSING/ChrS20.159.samp.0missing.vcf -grid 3541


