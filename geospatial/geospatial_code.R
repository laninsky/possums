# This code plots mitochondrial clades on a map

# 1. Loading necessary libraries
library(tidyverse)
library(ggmap)
library(ggrepel)

# 2. Setwd
setwd("Dropbox (Otago University)/PossumGenomeProject/PossumMitochondria/mtDNA_phylogenetic_analysis/")

# 3. Reading in data (tab delimited)
temp <- read_tsv("Codex_Possum_30Jun2022_mtDNA.txt")

# 4. Filtering to NZ points with lat/long
nz_geo <- temp %>% 
  filter(!(Location %in% c("lutruwita (Tasmania)", "Karta Pintingga (Kangaroo Island)","NSW"))) %>% 
  filter(!is.na(Longitude)) %>% 
  filter(!(-45.5<Latitude))

# 5. Making a box containing all sites along with some buffer
minlong <- min(nz_geo$Jitter_Long)-0.1
maxlong <- max(nz_geo$Jitter_Long)+0.1
minlat <- min(nz_geo$Jitter_Lat)-0.1
maxlat <- max(nz_geo$Jitter_Lat)+0.1

sbbox <- make_bbox(lon=c(minlong,maxlong), lat=c(minlat,maxlat),f=0)

# 6. Using the sbbox object to retrieve a map covering the sample sites
sq_map <- get_map(location = sbbox, maptype = "terrain-background", source = "stamen", crop=TRUE)

# 7. Creating map of sample locations by "mitochondrial clade"
# Pōhutukawa theme
# light blue, dark red,light green, dark green,
# c("#5FA1F7", "#9B1F1A","#83A552", "#3D4928")
# (c("BufferZone_005","Lawrence_004","NSW","TSW"))
ggmap(sq_map) +
  geom_point(data = nz_geo, mapping = aes(x = Jitter_Long, y = Jitter_Lat,fill = `RAxML run 1`), shape=21,color = "black",size=3) +
  scale_fill_manual(values=c("#5FA1F7", "#9B1F1A","#83A552", "#3D4928")) +
  coord_fixed(ratio=1) + xlab("Longitude") + ylab("Latitude") + 
  theme_bw(base_size = 20) +
  theme(legend.position="none",panel.border=element_rect(fill = NA)) +  
  theme(axis.title=element_text(size=28,face="bold"))

ggsave(filename="NZ_mitochondrial_site.pdf",plot = last_plot(),width=(maxlong-minlong)*10,height=(maxlat-minlat)*10,units="in")

# 8. Creating inserts for Lawrence and Woodside Glen to show clade distribution
# Using pies for these locations because of their more or less identical locations
Lawrence <- nz_geo %>% filter(Location=="Lawrence") %>% group_by(`RAxML run 1`) %>% summarise(count=n())

# Pōhutukawa theme
# light blue, dark red,light green, dark green,
# c("#5FA1F7", "#9B1F1A","#83A552", "#3D4928")
# (c("BufferZone_005","Lawrence_004","NSW","TSW"))

ggplot(Lawrence, aes(x = "", y = count, fill = `RAxML run 1`)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#5FA1F7", "#9B1F1A","#3D4928")) +
  theme_void() +
  theme(legend.position="none") +
  theme(plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

ggsave(filename="Lawrence.pdf",plot = last_plot(),width=2,height=2,units="in")

Woodside <- nz_geo %>% filter(Location %in% c("Woodside Glen","Woodside Glen?")) %>% group_by(`RAxML run 1`) %>% summarise(count=n())

# Pōhutukawa theme
# light blue, dark red,light green, dark green,
# c("#5FA1F7", "#9B1F1A","#83A552", "#3D4928")
# (c("BufferZone_005","Lawrence_004","NSW","TSW"))

ggplot(Woodside, aes(x = "", y = count, fill = `RAxML run 1`)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#83A552", "#3D4928")) +
  theme_void() +
  theme(legend.position="none") +
  theme(plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

ggsave(filename="Woodside.pdf",plot = last_plot(),width=2,height=2,units="in")

# 9. Creating inserts for Dunedin to show finer-scale distribution
Dunedin <- nz_geo %>% filter(!(Location %in% c("Lawrence", "Woodside Glen", "Woodside Glen?","Colony (Mother = Woodside Glen)"))) 

#Making a box containing all sites along with some buffer
minlong <- min(Dunedin$Jitter_Long)-0.1
maxlong <- max(Dunedin$Jitter_Long)+0.1
minlat <- min(Dunedin$Jitter_Lat)-0.1
maxlat <- max(Dunedin$Jitter_Lat)+0.1

sbbox <- make_bbox(lon=c(minlong,maxlong), lat=c(minlat,maxlat),f=0)

# Using the sbbox object to retrieve a map covering the sample sites
sq_map <- get_map(location = sbbox, maptype = "terrain-background", source = "stamen", crop=TRUE)

# Pōhutukawa theme
# light blue, dark red,light green, dark green,
# c("#5FA1F7", "#9B1F1A","#83A552", "#3D4928")
# (c("BufferZone_005","Lawrence_004","NSW","TSW"))

# Creating map of sample locations by "mitochondrial clade"
# blue, green, orange, red
ggmap(sq_map) +
  geom_point(data = Dunedin, mapping = aes(x = Jitter_Long, y = Jitter_Lat,fill = `RAxML run 1`), shape=21,color = "black",size=3) +
  scale_fill_manual(values=c("#5FA1F7","#83A552", "#3D4928")) +
  coord_fixed(ratio=1) + 
  theme(legend.position="none",panel.border=element_rect(fill = NA)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank() 
  )

ggsave(filename="Dunedin_jitter.pdf",plot = last_plot(),width=(maxlong-minlong)*10,height=(maxlat-minlat)*10,units="in")

ggmap(sq_map) +
  geom_point(data = Dunedin, mapping = aes(x = Longitude, y = Latitude,fill = `RAxML run 1`), shape=21,color = "black",size=3) +
  scale_fill_manual(values=c("#5FA1F7","#83A552", "#3D4928")) +
  coord_fixed(ratio=1) + 
  theme(legend.position="none",panel.border=element_rect(fill = NA)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank() 
  )

ggsave(filename="Dunedin_no_jitter.pdf",plot = last_plot(),width=(maxlong-minlong)*10,height=(maxlat-minlat)*10,units="in")


# Also do a summary pie for dunedin to make it consistent with others
Dunedin <- Dunedin %>% group_by(`RAxML run 1`) %>% summarise(count=n())

# Pōhutukawa theme
# light blue, dark red,light green, dark green,
# c("#5FA1F7", "#9B1F1A","#83A552", "#3D4928")
# (c("BufferZone_005","Lawrence_004","NSW","TSW"))

ggplot(Dunedin, aes(x = "", y = count, fill = `RAxML run 1`)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#5FA1F7","#83A552", "#3D4928")) +
  theme_void() +
  theme(legend.position="none") +
  theme(plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

ggsave(filename="Dunedin_pies.pdf",plot = last_plot(),width=2,height=2,units="in")

# 10. Sample sizes for figure
sum(Lawrence$count)
sum(Woodside$count)
sum(Dunedin$count)

# 11. Including Australian reference samples
temp <- temp %>% filter(!is.na(Longitude))

minlong <- 130
maxlong <- 180
minlat <- -30
maxlat <- -50

sbbox <- make_bbox(lon=c(minlong,maxlong), lat=c(minlat,maxlat),f=0)

# Using the sbbox object to retrieve a map covering the sample sites
sq_map <- get_map(location = sbbox, maptype = "terrain-background", source = "stamen", crop=TRUE)

# Pōhutukawa theme
#  taupe, light green, dark green
# c("#B19F8E","#83A552", "#3D4928")
# (c("KAN",NSW","TSW"))

ggmap(sq_map) +
  geom_point(data = temp %>% filter(Location %in% c("lutruwita (Tasmania)", "Karta Pintingga (Kangaroo Island)","NSW")), mapping = aes(x = Jitter_Long, y = Jitter_Lat,fill = `RAxML run 1`), shape=21,color = "black",size=3) +
  scale_fill_manual(values=c("#B19F8E","#83A552", "#3D4928")) +
  coord_fixed(ratio=1) + 
  theme(legend.position="none",panel.border=element_rect(fill = NA)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank() 
  )

ggsave(filename="Australasia.pdf",plot = last_plot(),width=(maxlong-minlong)/10,height=abs(maxlat-minlat)/10,units="in")
