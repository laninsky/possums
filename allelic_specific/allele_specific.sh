#!/bin/sh
echo "gVCF file: $1"
echo "RNA VCF file: $2"
echo "chromosome: $3"
echo "gene starting position: $4"
echo "gene ending position: $5"

echo "Subsetting gVCF file"
grep "^\#" "$1" >> $3.$1
grep -v "^\#" "$1" | grep "0\/1" | grep "$3" >> $3.$1

echo "Subsetting RNA VCF file"
grep "^\#" "$2" >> $3.$2
grep -v "^\#" "$2" | grep "$3" >> $3.$2
