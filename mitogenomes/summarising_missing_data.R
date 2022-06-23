# Loading libraries
library(tidyverse)
library(ggpubr)

# Directory with all fasta files in it
setwd("Dropbox (Otago University)/PossumGenomeProject/PossumMitochondria/mtDNA_phylogenetic_analysis/")

# For example:
# list.files(pattern="*.fasta")
#[1] "12SrRNA.fasta"                "16SrRNA.fasta"                "ATP6.fasta"                  
#[4] "ATP8.fasta"                   "COX1.fasta"                   "COX2.fasta"                  
#[7] "COX3.fasta"                   "CYTB.fasta"                   "DLOOP.fasta"                 
#[10] "ND1.fasta"                    "ND2.fasta"                    "ND3.fasta"                   
#[13] "ND4.fasta"                    "ND4L.fasta"                   "ND5.fasta"                   
#[16] "ND6.fasta"                    "origin_of_light_strand.fasta" "tRNA-Ala.fasta"              
#[19] "tRNA-Arg.fasta"               "tRNA-Asn.fasta"               "tRNA-Asp.fasta"              
#[22] "tRNA-Cys.fasta"               "tRNA-Gln.fasta"               "tRNA-Glu.fasta"              
#[25] "tRNA-Gly.fasta"               "tRNA-His.fasta"               "tRNA-Ile.fasta"              
#[28] "tRNA-Leu-1.fasta"             "tRNA-Leu-2.fasta"             "tRNA-Lys.fasta"              
#[31] "tRNA-Met.fasta"               "tRNA-Phe.fasta"               "tRNA-Pro.fasta"              
#[34] "tRNA-Ser-1.fasta"             "tRNA-Ser-2.fasta"             "tRNA-Thr.fasta"              
#[37] "tRNA-Trp.fasta"               "tRNA-Tyr.fasta"               "tRNA-Val.fasta" 

# Grouping files
tRNAs <- list.files(pattern="tRNA.*.fasta")
rRNAs <- list.files(pattern="rRNA.*.fasta")
other <- c("DLOOP.fasta","origin_of_light_strand.fasta")
coding <- list.files(pattern="fasta")[!(list.files(pattern="fasta") %in% c(tRNAs,rRNAs,other))]

# Summarizing missing data
initial_data <- NULL

for (i in tRNAs) {
  temp <- readLines(i)
  temp_names <- gsub(">","",temp[seq(1,length(temp),2)])
  temp_missing <- (nchar(temp[seq(2,length(temp),2)])*2-nchar(gsub("N","",temp[seq(2,length(temp),2)]))-nchar(gsub("-","",temp[seq(2,length(temp),2)],fixed=TRUE)))/nchar(temp[seq(2,length(temp),2)])*100 
  temp_out <- cbind(temp_names,temp_missing,i,"tRNA")
  initial_data <- rbind(initial_data,temp_out)
}

for (i in rRNAs) {
  temp <- readLines(i)
  temp_names <- gsub(">","",temp[seq(1,length(temp),2)])
  temp_missing <- (nchar(temp[seq(2,length(temp),2)])*2-nchar(gsub("N","",temp[seq(2,length(temp),2)]))-nchar(gsub("-","",temp[seq(2,length(temp),2)],fixed=TRUE)))/nchar(temp[seq(2,length(temp),2)])*100
  temp_out <- cbind(temp_names,temp_missing,i,"rRNA")
  initial_data <- rbind(initial_data,temp_out)
}

for (i in other) {
  temp <- readLines(i)
  temp_names <- gsub(">","",temp[seq(1,length(temp),2)])
  temp_missing <- (nchar(temp[seq(2,length(temp),2)])*2-nchar(gsub("N","",temp[seq(2,length(temp),2)]))-nchar(gsub("-","",temp[seq(2,length(temp),2)],fixed=TRUE)))/nchar(temp[seq(2,length(temp),2)])*100
  temp_out <- cbind(temp_names,temp_missing,i,"other")
  initial_data <- rbind(initial_data,temp_out)
}

for (i in coding) {
  temp <- readLines(i)
  temp_names <- gsub(">","",temp[seq(1,length(temp),2)])
  temp_missing <- (nchar(temp[seq(2,length(temp),2)])*2-nchar(gsub("N","",temp[seq(2,length(temp),2)]))-nchar(gsub("-","",temp[seq(2,length(temp),2)],fixed=TRUE)))/nchar(temp[seq(2,length(temp),2)])*100 
  temp_out <- cbind(temp_names,temp_missing,i,"coding")
  initial_data <- rbind(initial_data,temp_out)
}

# Takashi is completely missing for tRNA-Pro, so adding this missing data manually
initial_data <- rbind(initial_data,c("Takashi",68,"tRNA-Pro.fasta","tRNA"))

# Turning this data into a tibble
initial_data <- as_tibble(initial_data)
colnames(initial_data) <- c("Samples","Missing_data","Alignment","Region_type")
initial_data <- initial_data %>% mutate(Missing_data=as.numeric(Missing_data))

# Adding a column for "sample type"
initial_data <- initial_data %>% mutate(Sample_type=ifelse((Samples %in% c("Takashi", 
                                                          "NC_003039",
                                                          "RNA_posssum_ANU_liver",
                                                          "RNA_posssum_MWLR_liver",
                                                          "RNA_posssum_NSW_liver_1",
                                                          "RNA_posssum_NSW_liver_2",
                                                          "RNA_posssum_NSW_liver_3",
                                                          "RNA_posssum_NSW_liver_4",
                                                          "RNA_posssum_NSW_liver_5",
                                                          "RNA_posssum_NSW_liver_6",
                                                          "RNA_posssum_NSW_liver_7")), "Reference",
                                          ifelse(grepl("RNA",Samples),"RNA","LRPCR")))

# Plotting data for tRNA
ggplot(initial_data %>% filter(Region_type=="tRNA"),mapping=aes(x=Sample_type,y=Missing_data,colour=Sample_type)) + 
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(facets=vars(Alignment),scales="free")

ggsave("missing_tRNA.pdf",units = "in",width = 22,height=14)

# Plotting data for rRNA
ggplot(initial_data %>% filter(Region_type=="rRNA"),mapping=aes(x=Sample_type,y=Missing_data,colour=Sample_type)) + 
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(facets=vars(Alignment),scales="free")

ggsave("missing_rRNA.pdf",units = "in",width = 22,height=14)

# Plotting data for other
ggplot(initial_data %>% filter(Region_type=="other"),mapping=aes(x=Sample_type,y=Missing_data,colour=Sample_type)) + 
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(facets=vars(Alignment),scales="free")

ggsave("missing_other.pdf",units = "in",width = 22,height=14)

# Plotting data for coding
ggplot(initial_data %>% filter(Region_type=="coding"),mapping=aes(x=Sample_type,y=Missing_data,colour=Sample_type)) + 
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(facets=vars(Alignment),scales="free")

ggsave("missing_coding.pdf",units = "in",width = 22,height=14)

# Summarising missing data
mean_missing_data <- initial_data %>% group_by(Sample_type,Alignment,Region_type) %>% summarise(mean_missing=mean(Missing_data))

ggplot(mean_missing_data, mapping=aes(x=Sample_type,y=mean_missing,colour=Sample_type)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(facets=vars(Region_type),scales="free")

ggsave("missing__by_region.pdf",units = "in",width = 22,height=14)

mean_missing_data %>% group_by(Region_type,Sample_type) %>% summarise(max(mean_missing))
