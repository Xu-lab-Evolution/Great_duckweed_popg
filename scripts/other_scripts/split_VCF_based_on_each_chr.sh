#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 file.vcf(.gz) "
    exit 1
fi

VCF=$1

if [[ $VCF == *.gz ]]; then
    echo "VCF in gz format detected"
    VCFGZ=`basename $VCF`
else
    echo "VCF format detected"
    VCFGZ="${VCF##*/}.gz"
fi

dir=$(dirname $VCF)

if [ -f "$dir/$VCFGZ" ]; then
    echo "$dir/$VCFGZ exists, skip bgzipping..."
else
    echo "$dir/$VCFGZ dosn't exists, run bgzipping..."
    bgzip -c $VCF > $dir/$VCFGZ
fi


if [ -f "$dir/${VCFGZ}.tbi" ]; then
    echo "$dir/${VCFGZ}.tbi exists, skip tabix indexing..."
else
    echo "$dir/${VCFGZ}.tbi dosn't exists, run tabix indexing..."
    tabix -p vcf $dir/$VCFGZ
fi

tabix --list-chroms $dir/$VCFGZ > $dir/chromosomes.txt

while IFS= read -r line; do
  tabix -h $dir/$VCFGZ $line > $dir/$line.vcf;
done < $dir/chromosomes.txt
