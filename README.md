# possums
Leveraging hybrid zones to identify genes important to reproduction and survival in possums

### Mitogenomes
Received mitogenomes assembled to reference from RNAseq data from Oscar. Aligned these using the Geneious algorithm. Double-checked for "sane" protein-coding regions (e.g. no premature stop codons/frame-shift mutations etc). Checked entire alignment by eye and realigned indel regions for control region (only part of the alignment that looked a little more problematic).  

Following this, extracted the following partitions: concatenated rRNA genes, concatenated tRNA genes, concatenated protein-coding genes coded on the forward strand, ND6 (coded on the reverse strand), and control region/D-loop. Reverse-complemented ND6, and checked all protein-coding genes had nucleotides in multiples of threes, adding Ns if not (i.e. due to incomplete stop codons and/or removing of overlapping regions, as overlapping regions between any genes were excluded.
