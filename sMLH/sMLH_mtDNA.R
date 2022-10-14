# This code compares sMLH by sampling location and mitochondrial haplotype

# 1. Loading necessary libraries
library(tidyverse)
library(rstatix)

# 2. Setwd
setwd("Dropbox (Otago University)/PossumGenomeProject/PossumMitochondria/mtDNA_phylogenetic_analysis/")

# 3. Reading in data (tab delimited)
temp <- read_tsv("sMLH_vs_mtDNA.txt")
temp <- temp %>% select(sMLH,raxml_run1,Location,broad_location)

# 4. Filtering to just Otago and Australia
aus_otago <- temp %>% filter(broad_location %in% c("Otago","Australia"))

# Pōhutukawa theme
# light blue, taupe, light green, dark green,
# c("#5FA1F7", "#B19F8E", "#83A552", "#3D4928")
# (c("BufferZone_005", "KAN", "NSW","TSW"))

ggplot(aus_otago) + 
  geom_boxplot(mapping=aes(x=broad_location, y=sMLH,fill=raxml_run1)) +
  geom_jitter(mapping=aes(x=broad_location, y=sMLH,fill=raxml_run1),shape=21, position=position_jitter(0.2),size=1) +
  facet_wrap(~raxml_run1,nrow=1) +
  scale_fill_manual(values=c("#5FA1F7", "#B19F8E", "#83A552", "#3D4928")) +
  xlab("Location") +
  theme_bw(base_size = 8) +
  theme(legend.position="none",panel.border=element_rect(fill = NA)) +  
  theme(axis.title=element_text(size=14,face="bold")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave(filename="sMLH_mtDNA.pdf",plot = last_plot(),width=4,height=2,units="in")

#5. T-test within haplotypes and across all samples
t_test(data = aus_otago, 
       formula = sMLH ~ broad_location, 
       alternative = "two.sided")

t_test(data = aus_otago %>% filter(raxml_run1=="NSW"), 
       formula = sMLH ~ broad_location, 
       alternative = "two.sided")
