library(tidyverse)
library(cowplot)
library(readxl)

#### RNAseq data ####
load("data_clean/P337_BAL_data.RData")
load("data_clean/P337_BAL_module_data.RData")

#### Gene PCA ####
PCA <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  t() %>% 
  prcomp()

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1], digits=3)*100, "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2], digits=3)*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.BAL.abund.norm.voom$targets, by="libID")

PCA1 <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit, shape=visit),
             size=3) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, title="Genes") +
  coord_fixed(ratio=1) +
  scale_colour_manual(name = "SBP-Ag",
                      labels = c("Pre","Post"),
                      values = scales::hue_pal()(2)) +   
  scale_shape_manual(name = "SBP-Ag",
                     labels = c("Pre","Post"),
                     values = c(19,17)) +
  theme(legend.position = "bottom")
PCA1

#### Module PCA ####
PCA <- as.data.frame(mod.voom) %>% 
  column_to_rownames("module") %>% 
  t() %>% 
  prcomp()

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1], digits=3)*100, "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2], digits=3)*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.BAL.abund.norm.voom$targets, by="libID")

PCA2 <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit, shape=visit),
             size=3) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, title="Modules") +
  coord_fixed(ratio=1)  +
  theme(legend.position = "none")
  
#PCA2

#### Cytokine PCA ####
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("rowname","visit"), sep="_") %>% 
  mutate(name =paste(donorID, visit, sep="_")) %>% 
  select(-donorID, -visit) %>% 
  pivot_wider() %>% 
  column_to_rownames() %>% 
  t()

PCA <- plex %>% 
  prcomp()

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1], digits=3)*100, "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2], digits=3)*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("sampID") %>%
  # Select PCs
  dplyr::select(sampID, PC1:PC3) %>% 
  # Add  metadata
  separate(sampID, into=c("donorID","visit"), sep="_", remove = FALSE)

PCA4 <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit, shape=visit),
             size=3) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, title="Proteins") +
  coord_fixed(ratio=1)  +
  theme(legend.position = "none")

PCA4

#### Save ####
plot.R <- plot_grid(PCA2,PCA4, ncol=1,
                    labels = c("(B)", "(C)"), hjust=-2,
                    label_size = 13, label_fontface = "plain")
plot <- plot_grid(PCA1,plot.R, nrow=1,
          labels = c("(A)", ""), hjust=-1.2, vjust=1.8,
          label_size = 13, label_fontface = "plain")
#plot

ggsave("publication/FigS3.PCA.png", height=5, width=8)
ggsave("publication/FigS3.PCA.pdf", height=5, width=8)

