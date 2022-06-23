## Mitogenomes

### Alignment prep
Received mitogenomes assembled to reference from RNAseq data from Oscar and LR-PCR mitogenomes from Donna. Aligned these alignments together using the default Geneious algorithm. Checked entire alignment by eye and realigned indel regions. Conducted a quick and dirty UPGMA and NJ tree on the entire alignment to check duplicate samples of Sandy (e.g. different tissues from the reference individual) were clading together. 

Following this, I extracted annotated regions (e.g. tRNAs, rRNAs, origin of light strand, protein-coding regions, d-loop). Regions with no annotations, or overlapping annotations, were not extracted. I then double-checked for "sane" protein-coding regions (e.g. no premature stop codons/frame-shift mutations etc), and added extra Ns if protein-coding genes were not in multiples of three (for codon-partitioned models downstream). ND6 was reverse-complemented.

Following this, I exported these alignments to [evaluate levels of missingness](summarising_missing_data.R) for each region by data type (e.g. RNAseq vs LR-PCR). I found that the RNAseq samples had considerably more missing data than other data types for the majority of the tRNAs (all bar one), for both of the rRNAs, and for both of the "other" partitions (origin of light strand and d-loop). In contrast, only eight of the protein-coding genes had signficantly more missing data for the RNAseq samples (COX1, COX2, COX3, ND1, ND2, ND3, ND4L, ND5), and overall the average levels of missing data were lower (<1.24% for the coding gene with the most average missing data for any coding gene):

| Region type | Sample type | Maximum average missing data over all genes |
| ----------- | ----------- | ----------------------------------- |
| coding      | LRPCR       |              0.692                  |
| coding      | Reference   |              0.606                  |
| coding      | RNA         |              1.24                   |
| other       | LRPCR       |              3.67                   |
| other       | Reference   |              6.49                   |
| other       | RNA         |             33.6                    |
| rRNA        | LRPCR       |              0.0885                 |
| rRNA        | Reference   |              0.711                  |
| rRNA        | RNA         |              4.90                   |
| tRNA        | LRPCR       |             14.9                    |
| tRNA        | Reference   |             37.3                    |
| tRNA        | RNA         |             94.4                    |

Based on this, we decided to do phylogenetic reconstructions on the protein-coding regions only to limit the influence of how the mitogenome was generated for each sample. I next checked levels of missing data for each of the duplicated tissues of Sandy to [choose the sample with the least missing data](summarising_missing_data.R) across the protein-coding genes. We ended up retaining liver.

Following this, I concatenated the protein-coding genes. I then double-checked to make sure there were no discrepancies between our three maternally-related individuals, Sandy, Puku and Sheila (and there weren't!). I then partitioned out into codon 1, codon 2, and codon 3. All of these partitions were then concatenated, and the file was exported as [phylip](https://github.com/laninsky/possums/blob/main/mitogenomes/data/total_partitioned_alignment.phylip). The following file (`partion_file`) was generated based on alignment order/length for partitioning for the RAxML analysis:
```
DNA, PC_codon1 = 4086-7854
DNA, PC_codon2 = 7855-11623
DNA, PC_codon3 = 11624-15392
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
Following runs, downloaded RAxML_bipartitions file for each run and checked for topological convergence. After confirming this, run with the best likelihood (as presented in RAxML_info) was used as representative tree. No highly-supported discordant clades were found between the runs, and run 2 was found to have the highest likelihood. These outputs are available [here](https://github.com/laninsky/possums/tree/main/mitogenomes/output).

Although we initially [ran exabayes](old_analyses/exabayes.md) on previous exploratory analyses of the full mitogenome alignment, following the initial full-length runs (approximately 4 hours in length), ESS values associated with LnL and some parameters (particularly G-T transversions) were low for all partitions. This suggests that there is not enough variation across our partitions to implement GTR (the only model in ExaBayes). Because of this, for the current analyses where we've used even less of the full alignment (by only including the protein-coding genes), we exported our protein-coding partitions as nexus format so we could set up a BEAST analysis through BEAUTi.

## BEAST
Following the two BEAST runs, the logs were viewed in Tracer to determine appropriate burn-in. Following this, log and tree files had burn-in removed, and states thinned to leave approximately 20,000 states. Following this, TreeAnnotator was used to create a consensus tree for each run.
```
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.log -burnin 60 -resample 10000 -o beast_run2_thinned.log
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.trees -burnin 60 -resample 10000 -o beast_run2_thinned.trees
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/treeannotator -burnin 0 beast_run2_thinned.trees annotated_beast_run2.tre

/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.log -burnin 90 -resample 2000 -o beast_run1_thinned.log
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.trees -burnin 90 -resample 2000 -o beast_run1_thinned.trees
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/treeannotator -burnin 0 beast_run1_thinned.trees annotated_beast_run1.tre
