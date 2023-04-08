#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

gVCF_file <- args[1]
RNA_VCF_file <- args[2]
chromosome <- args[3]
start_pos <- as.numeric(args[4])
end_pos <- as.numeric(args[5])

# Lines for testing
gVCF_file <- "chr9.Sandy_S1_grouped.bam.vcf"
RNA_VCF_file <- "chr9.Sandy_liver.vcf"
start_pos <- 191485302
end_pos <- 191488654

# Library
library(tidyverse)
library(vcfR)
library(slider)

# Reading in the files
gVCF <- vcfR2tidy(read.vcfR(gVCF_file, verbose = FALSE))$gt
RNA_VCF <- vcfR2tidy(read.vcfR(RNA_VCF_file, verbose = FALSE))$gt 

# Getting minor allele frequency
gVCF_mod <- gVCF %>% filter(!is.na(gt_AD)) %>% mutate(gt_AD=gsub("(^.*)(,0$)","\\1",gt_AD)) %>% 
  rowwise() %>%  mutate(min_allele_freq=min(as.numeric(unlist(strsplit(gt_AD,","))))/sum(as.numeric(unlist(strsplit(gt_AD,",")))))

RNA_VCF_mod <- RNA_VCF %>% filter(!is.na(gt_AD)) %>% mutate(gt_AD=gsub("(^.*)(,0$)","\\1",gt_AD)) %>% 
  rowwise() %>%  mutate(min_allele_freq=min(as.numeric(unlist(strsplit(gt_AD,","))))/sum(as.numeric(unlist(strsplit(gt_AD,",")))))

# Sliding window for gVCF
gVCF_POS_mean_slide <- slide_dbl(gVCF_mod$POS, mean, .before = 10)
gVCF_min_allele_mean_slide <- slide_dbl(gVCF_mod$min_allele_freq, mean, .before = 10)
gVCF_min_allele_sd_slide <- slide_dbl(gVCF_mod$min_allele_freq, sd, .before = 10)
gVCF_slide <- as_tibble(cbind(gVCF_POS_mean_slide,gVCF_min_allele_mean_slide,gVCF_min_allele_sd_slide))
names(gVCF_slide) <- c("POS","mean_gVCF","sd_gVCF")
gVCF_slide <- gVCF_slide %>% mutate(min_gVCF=mean_gVCF-qt(0.975,df=9)*sd_gVCF/sqrt(10)) %>%
  mutate(max_gVCF=mean_gVCF+qt(0.975,df=9)*sd_gVCF/sqrt(10))

# Sliding window for RNA_VCF
RNA_VCF_POS_mean_slide <- slide_dbl(RNA_VCF_mod$POS, mean, .before = 10)
RNA_VCF_min_allele_mean_slide <- slide_dbl(RNA_VCF_mod$min_allele_freq, mean, .before = 10)
RNA_VCF_min_allele_sd_slide <- slide_dbl(RNA_VCF_mod$min_allele_freq, sd, .before = 10)
RNA_VCF_slide <- as_tibble(cbind(RNA_VCF_POS_mean_slide,RNA_VCF_min_allele_mean_slide,RNA_VCF_min_allele_sd_slide))
names(RNA_VCF_slide) <- c("POS","mean_RNA_VCF","sd_RNA_VCF")
RNA_VCF_slide <- RNA_VCF_slide %>% mutate(min_RNA_VCF=mean_RNA_VCF-qt(0.975,df=9)*sd_RNA_VCF/sqrt(10)) %>%
  mutate(max_RNA_VCF=mean_RNA_VCF+qt(0.975,df=9)*sd_RNA_VCF/sqrt(10))

# Plotting
ggplot() + 
  geom_ribbon(data=gVCF_slide,
              mapping=aes(x=POS,ymin=min_gVCF,ymax=max_gVCF),
              fill="#7d9fc2",alpha=0.5,colour="black") +
  geom_line(data=gVCF_slide,mapping=aes(x=POS,y=mean_gVCF),colour="#4d5f8e") +
  geom_point(data=RNA_VCF_mod,mapping=aes(x=POS,y=min_allele_freq),fill="#C582B2",size=5, shape=21) +
  geom_rect(aes(xmin=start_pos, ymin = 0, xmax = end_pos, ymax = 0.5),fill="#51806a",colour="black") +
  xlim((start_pos-100000),(end_pos+100000)) +
  theme_bw(base_size = 12) +
  xlab("Chromosomal position (bp)") +
  ylab("Minor allele frequency")

ggsave(paste(paste(chromosome, start_pos, end_pos, sep="_"),"pdf",sep="."))
