#!/bin/sh
echo "gVCF file: $1"
echo "RNA VCF file: $2"
echo "chromosome: $3"
echo "gene starting position: $4"
echo "gene ending position: $5"



echo "Subsetting RNA VCF file"
grep -v "^\#" "$2" | grep "$3" "$2" >> $3.$2
