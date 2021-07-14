# possums
Leveraging hybrid zones to identify genes important to reproduction and survival in possums

### Mitogenomes
Received mitogenomes assembled to reference from RNAseq data from Oscar. Aligned these using the default Geneious algorithm. Double-checked for "sane" protein-coding regions (e.g. no premature stop codons/frame-shift mutations etc). Checked entire alignment by eye and realigned indel regions for control region (only part of the alignment that looked a little more problematic). Conducted a quick and dirty UPGMA and NJ tree on the entire alignment. Based on these results, unclear whether there are two Tasmanian clades and one Mainland clade, or vice-versa. Due to confusion, decided to export partitioned mitogenome and try ML and Bayesian methods.  

Following this, extracted the following partitions: concatenated rRNA genes, concatenated tRNA genes, concatenated protein-coding genes coded on the forward strand, ND6 (coded on the reverse strand), and control region/D-loop. Reverse-complemented ND6, and checked all protein-coding genes had nucleotides in multiples of threes, adding Ns if not (i.e. due to incomplete stop codons and/or removing of overlapping regions, as overlapping regions between any genes were excluded. Protein-coding genes, including ND6, were then masked to specific codon positions. All partitions were then concatenated. The following file was generated for partitioning RAxML and ExaBayes analyses:
```
DNA, rRNA = 1-2510
DNA, tRNA = 2511-4085
DNA, PC_codon1 = 4086-7854
DNA, PC_codon2 = 7855-11623
DNA, PC_codon3 = 11624-15392
DNA, Dloop = 15393-17255
```
