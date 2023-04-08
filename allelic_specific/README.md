Need to have `allele_specific.sh` and `allele_specific.R` in the same folder as the `gVCF` and `RNA VCF file` given as arguments:
```
Usage is bash allele_specific.sh gVCF_file RNA_VCF_file chromosome gene_start_pos gene_end_pos
allele_specific.R and files should be in this directory along with allele_specific.sh
e.g.
bash allelic_specific.sh Sandy_S1_grouped.bam.vcf Sandy_liver.vcf chr9 191485302 191488654
```

R needs be in your path as well, as well as HarfBuzz
