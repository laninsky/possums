## Mitogenomes

### Alignment prep
Received mitogenomes assembled to reference from RNAseq data from Oscar and LR-PCR mitogenomes from Donna. Aligned these alignments together using the default Geneious algorithm. Checked entire alignment by eye and realigned indel regions. Conducted a quick and dirty UPGMA and NJ tree on the entire alignment to check duplicate samples of Sandy (e.g. different tissues from the reference individual) were clading together. 

Following this, I extracted annotated regions (e.g. tRNAs, rRNAs, origin of light strand, protein-coding regions, d-loop). Regions with no annotations, or overlapping annotations, were not extracted. I then double-checked for "sane" protein-coding regions (e.g. no premature stop codons/frame-shift mutations etc), and added extra Ns if protein-coding genes were not in multiples of three (for codon-partitioned models downstream). ND6 was reverse-complemented.

Following this, I exported these alignments to check levels of missing data per sample (to select for the representative tissue to use for Sandy), and to evaluate levels of missingness for each region by data type (e.g. RNAseq vs LR-PCR). We found that the RNAseq samples had considerably more missing data than other data types for the majority of the tRNAs (all bar one), for both of the rRNAs, and for both of the "other" partitions (origin of light strand and d-loop). In contrast, only eight of the protein-coding genes had signficantly more missing data for the RNAseq samples (COX1, COX2, COX3, ND1, ND2, ND3, ND4L, ND5), and overall the average levels of missing data were lower (<1.24% for the coding gene with the most average missing data for any coding gene):

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

Based on this, we decided to do phylogenetic reconstructions on the protein-coding regions only to limit the influence of how the mitogenome was generated for each sample.


# WHO WAS REMOVED BASED ON THIS, AND WERE ALL REGIONS USED?

Following this, I created six partitions: concatenated rRNA genes, concatenated tRNA genes, concatenated protein-coding gene codon 1 positions, codon 2 positions, codon 3 positions and control region/D-loop. All partitions were then concatenated, and file was exported as [phylip](https://github.com/laninsky/possums/blob/main/mitogenomes/data/total_partitioned_alignment.phylip). The following file (`partion_file`) was generated based on alignment order/length for partitioning RAxML and ExaBayes analyses:
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
#SBATCH -J exabayes_run1
#SBATCH --ntasks 1
#SBATCH -c 72
#SBATCH -t 48:00:00
#SBATCH --mem=105G
#SBATCH -D /nesi/nobackup/uoo03004/possums/mitogenomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1

module load ExaBayes/1.5.1-gimpi-2020a
yggdrasil -f total_partitioned_alignment.phylip -q exabayes_partition_file -s $RANDOM -n run1 -T 72 -M 0 -c config.nexus
```
```
#!/bin/bash -e

#SBATCH -A uoo03004
#SBATCH -J exabayes_run2
#SBATCH --ntasks 1
#SBATCH -c 72
#SBATCH -t 48:00:00
#SBATCH --mem=105G
#SBATCH -D /nesi/nobackup/uoo03004/possums/mitogenomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alana.alexander@otago.ac.nz
#SBATCH -N 1

module load ExaBayes/1.5.1-gimpi-2020a
yggdrasil -f total_partitioned_alignment.phylip -q exabayes_partition_file -s $RANDOM -n run2 -T 72 -M 0 -c config.nexus
```
Following the runs, checking parameter convergence:
```
module load ExaBayes/1.5.1-gimpi-2020a
postProcParam -n run1 -f ExaBayes_parameters*run1
postProcParam -n run2 -f ExaBayes_parameters*run2

# Checking overlap of parameters between runs in R
setwd("/Users/alanaalexander/Dropbox/Possums")
library(tidyverse)

# Reading in files 
data1 <- read_tsv("ExaBayes_parameterStatistics.run1")
data2 <- read_tsv("ExaBayes_parameterStatistics.run2")

# Creating a column to incorporate run
data1 <- data1 %>% mutate(run=1)
data2 <- data2 %>% mutate(run=2)

# Wrapping these in to one file
data <- rbind(data1,data2)

data
## A tibble: 138 x 11
#   paramName       mean       sd    perc5   perc25   median    perc75     per95    ESS  psrf   run
#   <chr>          <dbl>    <dbl>    <dbl>    <dbl>    <dbl>     <dbl>     <dbl>  <dbl> <dbl> <dbl>
# 1 alpha{5}     0.0666   0.0506   0.0544   0.0273   0.0544    0.0896    0.156   4295.   1.00     1
# 2 alpha{3}    97.5     58.7     96.5     46.2     96.5     149.      189.      5004.   1.00     1
# 3 alpha{2}    94.8     59.8     93.7     42.3     93.7     147.      190.      4816.   1.00     1
# 4 alpha{1}    98.9     58.3     99.9     47.6     99.9     148.      191.      5064.   1.00     1
# 5 r{5}(G<->T)  0.00465  0.00316  0.00366  0.00209  0.00366   0.00659   0.0107    20.3  1.27     1
# 6 r{5}(A<->C)  0.0140   0.0104   0.0116   0.00626  0.0116    0.0192    0.0342  3275.   1.00     1
# 7 r{4}(G<->T)  0.00413  0.00267  0.00337  0.00194  0.00337   0.00584   0.00924   19.6  1.28     1
# 8 r{5}(C<->T)  0.340    0.0559   0.338    0.301    0.338     0.379     0.435   2166.   1.00     1
# 9 r{4}(C<->T)  0.300    0.0355   0.299    0.275    0.299     0.324     0.360   1604.   1.00     1
#10 r{4}(C<->G)  0.0462   0.0276   0.0413   0.0256   0.0413    0.0614    0.0995  2464.   1.00     1
## … with 128 more rows

for (i in unique(data$paramName)) {
  tempdata <- data[which(data$paramName==i),]
  ggplot(tempdata) +
    geom_boxplot(mapping=aes(x=as.factor(run), ymin=perc5, lower=perc25, middle=median, upper=perc75, max=per95),  stat = "identity") + ggtitle(i) + theme_bw() + theme(
      plot.title = element_text(hjust=0.5))
  ggsave(paste(i,".pdf"))
}

data %>% filter(ESS<=200) %>% select(paramName,ESS,run) %>% arrange(paramName)
```
Following the initial full-length runs (approximately 4 hours in length), ESS values associated with LnL and some parameters (particularly G-T transversions) were low for all partitions. This suggests that there is not enough variation across our partitions to implement GTR (the only model in ExaBayes). Because of this, we went back to our partition alignments in Geneious, and exported them as nexus format so we could set up a BEAST analysis through BEAUTi.

## BEAST
Following the two BEAST runs, the logs were viewed in Tracer to determine appropriate burn-in. Following this, log and tree files had burn-in removed, and states thinned to leave approximately 20,000 states. Following this, TreeAnnotator was used to create a consensus tree for each run.
```
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.log -burnin 60 -resample 10000 -o beast_run2_thinned.log
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.trees -burnin 60 -resample 10000 -o beast_run2_thinned.trees
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/treeannotator -burnin 0 beast_run2_thinned.trees annotated_beast_run2.tre

/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.log -burnin 90 -resample 2000 -o beast_run1_thinned.log
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/logcombiner -log PC_codon3.trees -burnin 90 -resample 2000 -o beast_run1_thinned.trees
/opt/nesi/CS400_centos7_bdw/BEAST/2.6.3/bin/treeannotator -burnin 0 beast_run1_thinned.trees annotated_beast_run1.tre
