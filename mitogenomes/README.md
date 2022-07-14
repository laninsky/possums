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

Following this, I concatenated the protein-coding genes. I then double-checked to make sure there were no discrepancies between our three maternally-related individuals, Sandy, Puku and Sheila (and there weren't!). I then partitioned out into codon 1, codon 2, and codon 3. All of these partitions were then concatenated, and the file was exported as [phylip](data/concatenated_codon_partitioned_protein_coding_genes.phy). The following file (`partion_file`) was generated based on alignment order/length for partitioning for the RAxML analysis:
```
DNA, PC_codon1 = 1-3761
DNA, PC_codon2 = 3762-7522
DNA, PC_codon3 = 7523-11283
```

## RAxML
Ran two separate RAxML runs to check convergence. Used GTRCAT on advice of RaxML author. Dialed in resource use for run 2 based on the `nn_seff` results from run 1.
```
#!/bin/bash -e

#SBATCH -A uoo03398 
#SBATCH -J raxml_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 24:00:00
#SBATCH --mem=1G
#SBATCH -D /nesi/nobackup/uoo03398/possum_mitogenomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load RAxML/8.2.12-gimkl-2020a
raxmlHPC-PTHREADS-SSE3 -s concatenated_codon_partitioned_protein_coding_genes.phy -q partition_file -n run1 -m GTRCAT -f a -N 100 -x $RANDOM -p $RANDOM -T 12
```
```
#!/bin/bash -e

#SBATCH -A uoo03398 
#SBATCH -J raxml_run2
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 2:00:00
#SBATCH --mem=200M
#SBATCH -D /nesi/nobackup/uoo03398/possum_mitogenomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load RAxML/8.2.12-gimkl-2020a
raxmlHPC-PTHREADS-SSE3 -s concatenated_codon_partitioned_protein_coding_genes.phy -q partition_file -n run2 -m GTRCAT -f a -N 100 -x $RANDOM -p $RANDOM -T 12
```
Following runs, downloaded RAxML_bipartitions file for each run and checked for topological convergence. After confirming this, run with the best likelihood (as presented in RAxML_info) was used as representative tree. No highly-supported discordant clades were found between the runs, and run 1 was found to have the highest likelihood. These outputs are available [here](https://github.com/laninsky/possums/tree/main/mitogenomes/outputs).

Although we [ran exabayes](old_analyses/exabayes.md) on previous exploratory analyses of the full mitogenome alignment, following the initial full-length runs (approximately 4 hours in length), ESS values associated with LnL and some parameters (particularly G-T transversions) were low for all partitions. This suggests that there is not enough variation across our partitions to implement GTR (the only model in ExaBayes). Because of this, for the current analyses where we've used even less of the full alignment (by only including the protein-coding genes), we exported our protein-coding partitions as nexus format so we could set up a Bayesian BEAST analysis through BEAUTi.

## BEAST
In our previous explatory analyses of the full mitogenome alignment we progressively simplified the site model in rensponse to poor ESS values, eventually settling on a gamma site model (4 categories, initial shape 1.0, estimated) and TN93 with Kappa1 and Kappa2 set at 2.0, but estimated. Site frequencies were also estimated. We did not estimate proportion invariant due to the interactions that can set up between that and a gamma site model. We implemented a Strict Clock, as we did not expect marked rate variation given all individuals were from the same species. We used a Coalescent Bayesian Skyline as our tree model due to uncertainty about population demography. Tree was linked across all partitions, but site and clock was allowed to vary. We ran two runs of 500,000,000, sampling every 10,000 states. We used these same parameters for our protein-coding gene only analyses (partitioned by codon).
```
#!/bin/bash -e

#SBATCH -A uoo03398 
#SBATCH -J beast_run1
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 72:00:00
#SBATCH --mem=8G
#SBATCH -D /nesi/nobackup/uoo03398/possum_mitogenomes/beast/run1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load BEAST/2.6.6
beast -threads 12 beast_possum_mitogenome_Jun2022.xml

#!/bin/bash -e

#SBATCH -A uoo03398 
#SBATCH -J beast_run2
#SBATCH --ntasks 1
#SBATCH -c 12
#SBATCH -t 72:00:00
#SBATCH --mem=8G
#SBATCH -D /nesi/nobackup/uoo03398/possum_mitogenomes/beast/run2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1
#SBATCH --hint=nomultithread

module load BEAST/2.6.6
beast -threads 12 beast_possum_mitogenome_Jun2022.xml

```
The BEAST runs timed out between 250,000,000 and 290,000,000, but after inspecting the logs in Tracer, these chains were more than long enough to reach decent (>200) ESS values. So I trimmed all the output files to 250,000,000 states for both runs. 
```
head -n 25135 run1.log > run1_trimmed.log
head -n 25135 run2.log > run2_trimmed.log
head -n 25661 run1.trees > run1_trimmed.trees
head -n 25661 run2.trees > run2_trimmed.trees
```
The log files were then checked again in Tracer to ensure adequate ESS values (all >200 except run1 bPopSizes.1 = 162.3, however converged on estimate of run2 so proceeded) and convergence between runs. Based on Tracer, the "usual" burn-in of 10% looked appropriate. Following this, TreeAnnotator was used to create a consensus tree for each run (Maximum clade credibility tree, Common Ancestor heights), removing the first 10% as burn-in.

The tree topology was then compared for each run. Following confirmation there were no highly supported differences, the two runs were combined, checked for adequate ESS values, and then a combined consensus tree was created and used as the representative Bayesian tree for this analysis.
