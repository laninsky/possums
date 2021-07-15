## Mitogenomes

### Alignment prep
Received mitogenomes assembled to reference from RNAseq data from Oscar. Aligned these using the default Geneious algorithm. Double-checked for "sane" protein-coding regions (e.g. no premature stop codons/frame-shift mutations etc). Checked entire alignment by eye and realigned indel regions for control region (only part of the alignment that looked a little more problematic). Conducted a quick and dirty UPGMA and NJ tree on the entire alignment. Based on these results, unclear whether there are two Tasmanian clades and one Mainland clade, or vice-versa. Due to confusion, decided to export partitioned mitogenome and try ML and Bayesian methods.  

Following this, extracted the following partitions: concatenated rRNA genes, concatenated tRNA genes, concatenated protein-coding genes coded on the forward strand, ND6 (coded on the reverse strand), and control region/D-loop. Reverse-complemented ND6, and checked all protein-coding genes had nucleotides in multiples of threes, adding Ns if not (i.e. due to incomplete stop codons and/or removing of overlapping regions, as overlapping regions between any genes were excluded. Protein-coding genes, including ND6, were then masked to specific codon positions. All partitions were then concatenated, and file was exported as [phylip](https://github.com/laninsky/possums/blob/main/mitogenomes/data/total_partitioned_alignment.phylip). The following file (`partion_file`) was generated based on alignment order/length for partitioning RAxML and ExaBayes analyses:
```
DNA, rRNA = 1-2510
DNA, tRNA = 2511-4085
DNA, PC_codon1 = 4086-7854
DNA, PC_codon2 = 7855-11623
DNA, PC_codon3 = 11624-15392
DNA, Dloop = 15393-17255
```

## RAxML
Ran two separate RAxML runs to check convergence. used GTRCAT on advice of RaxML author (initially ran run1 with `#SBATCH --qos=debug` and time limited to 15:00, with the intention of increasing the time once I confirmed that the runs were working. However, runs completed in ~2 minutes).
```
#!/bin/bash -e

#SBATCH -A uoo03004 
#SBATCH -J raxml_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 15:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo03004/possums/mitogenomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --qos=debug

module load RAxML/8.2.12-gimkl-2020a
raxmlHPC-PTHREADS-SSE3 -s total_partitioned_alignment.phylip -q partition_file -n run1 -m GTRCAT -f a -N 100 -x $RANDOM -p $RANDOM -T 12
```
```
#!/bin/bash -e

#SBATCH -A uoo03004 
#SBATCH -J raxml_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 15:00
#SBATCH --mem=20G
#SBATCH -D /nesi/nobackup/uoo03004/possums/mitogenomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread
#SBATCH --qos=debug

module load RAxML/8.2.12-gimkl-2020a
raxmlHPC-PTHREADS-SSE3 -s total_partitioned_alignment.phylip -q partition_file -n run1 -m GTRCAT -f a -N 100 -x $RANDOM -p $RANDOM -T 12
```
Following runs, downloaded RAxML_bipartitions file for each run and checked for topological convergence. After confirming this, run with the best likelihood (as presented in RAxML_info) used as representative tree, following Alexander and Short (XXXX). No highly-supported discordant clades were found between the runs, and run 2 was found to have the highest likelihood. These outputs are available [here](https://github.com/laninsky/possums/tree/main/mitogenomes/output).

## ExaBayes
To run ExaBayes, created the config.nexus file exabayes needs (contents below), following the methods of https://github.com/laninsky/beetles/tree/master/hydrophiloidea_analyses:
```
begin run; 
   numRuns 4
   numCoupledChains 3
end;
```
Also slightly tweaked `partition_file` to `exabayes_partition_file` by removing the spaces before and after the = sign:
```
DNA, rRNA=1-2510
DNA, tRNA=2511-4085
DNA, PC_codon1=4086-7854
DNA, PC_codon2=7855-11623
DNA, PC_codon3=11624-15392
DNA, Dloop=15393-17255
```
Ran two separate runs to check convergence (initially ran run1 with `#SBATCH --qos=debug` and time limited to 15:00, with the intention of increasing the time once I confirmed that the runs were working).
```
#!/bin/bash -e

#SBATCH -A uoo03004
#SBATCH -J 50perc_run1
#SBATCH --ntasks 1
#SBATCH -c 36
#SBATCH -t 15:00
#SBATCH --mem=105G
#SBATCH -D /nesi/nobackup/uoo00105/beetles/50perc_raxml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --qos=debug

module load ExaBayes/1.5.1-gimpi-2020a
exabayes -f total_partitioned_alignment.phylip -q exabayes_partition_file -s $RANDOM -n run1 -T 72 -M 0 -c config.nexus
```
