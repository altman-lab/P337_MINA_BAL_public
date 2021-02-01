library(tidyverse)
library(cowplot)

#### Data ####
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
  theme(legend.position = "none")
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
  coord_fixed(ratio=1) +
  scale_color_discrete(labels = c("Pre", "Post"),
                       name = c("Allergen challenge")) +
  theme(legend.position = "top")
  
#PCA2

#### Save ####

plot <- plot_grid(PCA1,PCA2, ncol=1,
          labels = c("(A) Genes", "(B) Modules"), hjust = 0)
plot

ggsave("publication/FigX.PCA.png", height=7, width=4.5)
