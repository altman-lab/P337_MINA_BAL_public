library(tidyverse)
library(readxl)
library(Hmisc)
library(patchwork)

#### Data ####
load("data_clean/P337_BAL_data.RData")

### cell % ####
##### Stats #####
cell.test <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, ends_with(".pct")) %>% 
  pivot_longer(-c(donorID, visit)) %>% 
  arrange(name, donorID, visit) %>% 
  group_by(name, visit) %>% 
  summarise(value = list(value)) %>% 
  pivot_wider(names_from = visit) %>% 
  group_by(name) %>% 
  mutate(p = t.test(unlist(V4), unlist(V5))$p.value,
         t = t.test(unlist(V4), unlist(V5))$statistic) %>% 
  ungroup() %>% 
  mutate(FDR = p.adjust(p))

##### plot #####
p1 <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, ends_with(".pct")) %>% 
  pivot_longer(-c(donorID, visit)) %>% 
  group_by(name, donorID) %>% 
  mutate(diff = value[visit=="V5"]-value[visit=="V4"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) %>%
  ungroup() %>% 
  left_join(cell.test) %>% 
  mutate(cell = gsub(".pct","",name),
         cell = recode(cell, 
                       "MONO"="monocyte",
                       "EOS"="eosinophil", 
                       "LYM"="lymphocyte",
                       "Epi"="epithelial",
                       "NEUT"="neutrophil"),
         FDR.format = ifelse(FDR > 0.01, round(FDR, digits=3),
                             formatC(FDR, digits = 2, format="e")),
         cell = paste(cell, paste("FDR =", FDR.format), sep="\n"),
         visit = recode_factor(factor(visit), "V4"="Pre", "V5"="Post")) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic(base_size = 8) +
  labs(x="Segmented bronchial provocation with allergen (SBP-Ag)",
       y="Cell %", color="Post - Pre\nchange") +
  facet_grid(~cell, scales="free_x") +
  scale_color_manual(values=c("down"="#74add1","up"="#FF8679"))
# p1

#### fMRI ####
#Time point 1 = visit 4
neuro <- read_excel(sheet="T1",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(visit="V4")

#Time point 2 = visit 5
neuro <- read_excel(sheet="T2",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(visit="V5") %>% 
  #Combine with other visit
  full_join(neuro)  %>% 
  #Format idnum to match RNAseq data
  mutate(idnum = paste("MA",idnum, sep="")) %>% 
  rename(donorID=idnum) %>% 
  #remove non fMRI vars
  filter(neuro != "BDI" & neuro != "LSI") %>% 
  #remove fxnal metrics
  filter(!grepl("3way",neuro) & !grepl("LPR", neuro)) %>% 
  #Remove failed sample
  filter(donorID != "MA1012") %>% 
  #Clean names
  mutate(neuro = gsub("_"," ",neuro)) %>%  
  mutate(visit = recode_factor(factor(visit), "V4"="Pre","V5"="Post")) %>% 
  mutate(neuro = gsub("amgy","amyg",neuro))

# Calculate delta
neuro.delta <- neuro %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=Post-Pre) %>% 
  #shorten donor ID
  mutate(donorID = gsub("^MA10","", donorID))

#### clustering ####
#Cluster order
neuro.delta.wide <- neuro.delta  %>%
  #wide format
  dplyr::select(-Pre,-Post) %>%
  pivot_wider(names_from = neuro, values_from = delta) %>%
  column_to_rownames("donorID")

ord.donor <- neuro.delta.wide %>% 
  dist(., method = "euclidean") %>% 
  hclust(., method="ward.D")

ord.neuro <- t(neuro.delta.wide) %>% 
  dist(., method = "euclidean") %>% 
  hclust(., method="ward.D")

#### Heatmap delta ####
p2 <- neuro.delta %>% 
  mutate(donorID.ord = factor(donorID, 
                              levels=rownames(neuro.delta.wide)[ord.donor$order]),
         neuro.ord = factor(neuro, 
                            levels=colnames(neuro.delta.wide)[ord.neuro$order])) %>%
  
  ggplot(aes(x=donorID.ord, y=neuro.ord, fill=delta)) +
  geom_tile() +
  labs(x="Donor", fill="Post - Pre BOLD score") +
  # coord_flip() +
  coord_fixed()+
  theme(text = element_text(size=8),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top", legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5,
                               title.position = "top", title.hjust = 0.5)) +
  scale_fill_gradient2(low="darkblue",mid="white",high="darkred",
                       limits=c(-105,105))
# p2

#### Cell data ####
cell.delta <- dat.BAL.abund.norm.voom$targets %>% 
  select(visit, donorID, ends_with(".pct")) %>% 
  pivot_longer(EOS.pct:Epi.pct) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=V5-V4) %>% 
  select(-V4, -V5) %>% 
  pivot_wider(values_from = delta) %>% 
  arrange(donorID) %>% 
  #shorten donor ID
  mutate(donorID = gsub("^MA10","", donorID))

#### combine data ####
X <- neuro.delta %>% 
  select(donorID, neuro, delta) %>% 
  pivot_wider(names_from = neuro, values_from = delta) %>% 
  full_join(cell.delta) %>%
  rename_all(~gsub(".pct", " %",.)) %>%
  # rename_all(~gsub("_"," ",.)) %>%
  column_to_rownames("donorID") %>%
  as.matrix()

#### Correlation ####
cor.result <- Hmisc::rcorr(X, type = "pearson")

R <- as.data.frame(cor.result$r) %>% 
  rownames_to_column() %>% 
  filter(grepl(" %",rowname)) %>% 
  column_to_rownames() %>% 
  select(!ends_with(" %")) %>% 
  as.matrix()

P <- as.data.frame(cor.result$P) %>% 
  rownames_to_column() %>% 
  filter(grepl(" %",rowname)) %>% 
  column_to_rownames() %>% 
  select(!ends_with(" %")) %>% 
  as.matrix()

#### clustering ####
#Cluster order
ord.cell <- R %>% 
  dist(., method = "euclidean") %>% 
  hclust(., method="ward.D")

#Use same neuro order as first heatmap

#### combine all pearson data ####
cor.arrange <- as.data.frame(cor.result$r) %>% 
  rownames_to_column("X_variable") %>% 
  pivot_longer(-X_variable, values_to = "r")
cor.arrange <- as.data.frame(cor.result$P) %>% 
  rownames_to_column("X_variable") %>% 
  pivot_longer(-X_variable, values_to = "p") %>% 
  full_join(cor.arrange) %>% 
  filter(grepl(" %",X_variable) & !grepl(" %",name)) %>% 
  #order by cluster
  mutate(cell.ord = factor(X_variable, 
                           levels=rownames(R)[ord.cell$order]),
         neuro.ord = factor(name, 
                            levels=colnames(neuro.delta.wide)[ord.neuro$order])) %>%
  #signif label
  mutate(label = ifelse(p < 0.05, "*",""),
         FDR = p.adjust(p),
         label2 = ifelse(FDR < 0.05, "*","")) %>% 
  #cell labels
  mutate(cell.lab = recode(X_variable,
                           "EOS %"="eosinophil",
                           "NEUT %"="neutrophil",
                           "LYM %"="lymphocyte",
                           "MONO %"="monocyte",
                           "Epi %"="epithelial"))

#### Heatmap cell % pearson ####
p3 <- cor.arrange %>% 
  
  ggplot(aes(x=reorder(cell.lab, cell.ord), y=neuro.ord, fill=r)) +
  geom_tile() +
  geom_text(aes(label = label)) +
  labs(x="Cell %", fill="Pearson (r)") +
  # coord_flip() +
  coord_fixed()+
  theme(text = element_text(size=8),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "top", legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5,
                               title.position = "top", title.hjust = 0.5)) +
  scale_fill_gradient2(low="darkblue",mid="white",high="darkred",
                       limits=c(-1,1))
# p3


#### Save ####

plot <- p1 / (p2 + p3) + 
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(heights = c(1,1.4))

# plot

#### Save ####
ggsave("publication/FigS1.cells.fMRI.png", plot, height=4, width=6)
ggsave("publication/FigS1.cells.fMRI.pdf", plot, height=4, width=6)
