## sMLH against mitochondrial haplotype

Oscar ran sMLH on nuclear SNPs of samples characterised using RNAseq in our dataset. I then combined this with the mitochondrial haplotype for those samples to see if there were any patterns by haplotype/location and by possums sampled in the introduced vs native range.

The necessary columns for [these analyses](sMLH_mtDNA.R) are as follows:
- sMLH # The sMLH for each individual as calculated by Oscar (currently plotted the 95% complete SNP sMLH calculations)
- raxml_run1 # The mitochondrial clade for each sample
- Location # The sampling locality for each sample, if known
- broad_location # Because our focus is on Otago individuals, OTAGO, NZ, or AUSTRALIA

e.g. 
```
# A tibble: 91 × 4
    sMLH raxml_run1      Location       broad_location
   <dbl> <chr>           <chr>          <chr>         
 1 0.940 Blue (NSW)      NA             OTAGO         
 2 1.11  Green (Takashi) NA             OTAGO         
 3 0.924 Green (Takashi) NA             OTAGO         
 4 0.774 Blue (NSW)      NA             OTAGO         
 5 1.02  Red (mixed)     South Dale Rd? OTAGO         
 6 0.826 Blue (NSW)      NA             OTAGO         
 7 0.939 Green (Takashi) South Dale Rd  OTAGO         
 8 0.673 Green (Takashi) South Dale Rd  OTAGO         
 9 0.948 Blue (NSW)      Woodside Glen  OTAGO         
10 1.08  Green (Takashi) South Dale Rd  OTAGO         
# … with 81 more rows
```
