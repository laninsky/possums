#!/bin/sh

echo "Usage is bash allele_specific.sh gVCF_file RNA_VCF_file chromosome gene_start_pos gene_end_pos"
echo "allele_specific.R and files should be in this directory with allele_specific.sh"
echo "e.g."
echo " bash allelic_specific.sh Sandy_S1_grouped.bam.vcf Sandy_liver.vcf chr9 191485302 191488654"

echo "gVCF file: $1"
echo "RNA VCF file: $2"
echo "chromosome: $3"
echo "gene starting position: $4"
echo "gene ending position: $5"

if [ -f "$3.$1" ]; then
  echo "$3.$1 exists. Proceeding to next step"
else
  echo "Subsetting gVCF file"
  grep "^\#" "$1" >> $3.$1
  grep -v "^\#" "$1" | grep "0\/1" | grep "$3" >> $3.$1
fi

if [ -f "$3.$2" ]; then
  echo "$3.$2 exists. Proceeding to next step"
else
  echo "Subsetting RNA VCF file"
  grep "^\#" "$2" >> $3.$2
  grep -v "^\#" "$2" | grep "$3" >> $3.$2
fi
