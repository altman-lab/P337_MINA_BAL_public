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
  geom_point(aes(color=visit),
             size=3) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label) +
  coord_fixed(ratio=1) +
  scale_color_discrete(labels = c("Pre", "Post"),
                       name = c("Allergen challenge")) +
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
  geom_point(aes(color=visit),
             size=3) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label) +
  coord_fixed(ratio=1)  +
  theme(legend.position = "none")
  
#PCA2

#### fMRI data ####
#Time point 1 = visit 4
neuro <- read_excel(sheet="T1",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(sampID=paste("MA", idnum,"_V4",sep=""))

#Time point 2 = visit 5
neuro <- read_excel(sheet="T2",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(sampID=paste("MA", idnum, "_V5",sep=""))  %>% 
  #Combine with other visit
  full_join(neuro) %>% 
  #wide format
  pivot_wider(names_from = neuro) %>% 
  select(-BDI, -LSI, -idnum) %>% 
  column_to_rownames("sampID")

#### fMRI PCA ####
PCA <- neuro %>% 
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

PCA3 <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit),
             size=3) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label) +
  coord_fixed(ratio=1)  +
  theme(legend.position = "none")

#PCA3

#### Save ####
plot.R <- plot_grid(PCA2,PCA3, ncol=1,
                    labels = c("(B) Modules", "(C) fMRI"), hjust=0)
plot <- plot_grid(PCA1,plot.R, nrow=1,
          labels = c("(A) Genes", ""), hjust=0)
plot

ggsave("publication/FigX.PCA.png", height=6, width=8)
