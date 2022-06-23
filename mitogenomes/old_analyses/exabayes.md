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
