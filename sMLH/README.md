## sMLH against mitochondrial haplotype

Oscar ran sMLH on nuclear SNPs of samples characterised using RNAseq in our dataset. I then combined this with the mitochondrial haplotype for those samples to see if there were any patterns by haplotype/location and by possums sampled in the introduced vs native range.

The necessary columns for [these analyses](sMLH_mtDNA.R) are as follows:
- sMLH # The sMLH for each individual as calculated by Oscar (currently plotted the 95% complete SNP sMLH calculations)
- raxml_run1 # The mitochondrial clade for each sample
- Location # The sampling locality for each sample, if known
- broad_location # Because our focus is on Otago individuals, OTAGO, NZ, or AUSTRALIA
