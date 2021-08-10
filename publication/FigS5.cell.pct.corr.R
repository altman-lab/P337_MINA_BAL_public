library(tidyverse)
library(corrplot)
library(readxl)
library(Hmisc)

#### Cell data ####
load("data_clean/P337_BAL_data.RData")

cell.delta <- dat.BAL.abund.norm.voom$targets %>% 
  select(visit, donorID, EOS.pct, PMN.pct) %>% 
  pivot_longer(EOS.pct:PMN.pct) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=V5-V4) %>% 
  select(-V4, -V5) %>% 
  pivot_wider(values_from = delta)
  
#### fMRI data ####
#Time point 1 = visit 4
neuro <- read_excel(sheet="T1",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add visit variable
  mutate(visit="V4")

#Time point 2 = visit 5
neuro <- read_excel(sheet="T2",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add visit variable
  mutate(visit="V5") %>% 
  #Combine with other visit
  full_join(neuro) %>% 
  #Format idnum to match RNAseq data
  mutate(idnum = paste("MA",idnum, sep="")) %>% 
  rename(donorID=idnum)

neuro.delta <- neuro %>% 
  filter(donorID != "MA1012" & neuro != "LSI" & neuro != "BDI") %>% 
  #Calculate delta
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  #remove fxnal metrics
  filter(!grepl("3way",neuro) & !grepl("LPR", neuro)) %>% 
  #fix name
  mutate(neuro = gsub("amgy","amyg",neuro)) %>% 
  #wide format
  dplyr::select(-V4,-V5) %>% 
  pivot_wider(names_from = neuro, values_from = delta) 

##### Module data ####
attach("data_clean/P337_BAL_module_data.RData")

mod.delta <- mod.voom %>% 
  pivot_longer(-module, names_to = "libID") %>% 
  left_join(select(dat.BAL.abund.norm.voom$targets, libID, visit, donorID)) %>% 
  select(-libID) %>% 
  filter(grepl("EOS.pct_02", module)) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(EOS02=V5-V4) %>% 
  select(-V4, -V5, -module)

##### protein data ####
#SPLS-selected
attach("results/PLS/SPLS.RData")

cyto <- as.data.frame(result.spls$loadings$X) %>% 
  filter(comp1 !=0 | comp2 != 0)

plex.delta<- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("name","visit")) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=V5-V4) %>% 
  select(-V4, -V5) %>% 
  filter(name %in% rownames(cyto)) %>% 
  pivot_wider(values_from = delta) %>% 
  rename(donorID=ptID)

#### combine data ###
X <- full_join(cell.delta, mod.delta) %>% 
  full_join(plex.delta) %>% 
  full_join(neuro.delta) %>% 
  rename(`EOS02 module`=EOS02, Eosinophil=EOS.pct, PMN=PMN.pct) %>% 
  rename_all(~gsub("_"," ",.)) %>% 
  #order Y
  select(donorID, `Eosinophil`, `PMN`,
         `FLT3L`,`IL17A`,`EOS02 module`,`GCSF`,`IL10`,
         `PDGFAA`,`CCL27`,`IFNA2`,`GMCSF`,`CCL21`,
         `TGFA`,`CCL26`,`IL23`,`CCL7`,`CXCL1`,`IL16`,
         `L vA insula`,`L HO amyg`,`L dA insula`,`R HO amyg`,
         `pACC`,`dACC`,`R dA insula`,`R vA insula`) %>% 
  column_to_rownames("donorID") %>% 
  as.matrix()
# 
# Y <- cell.delta %>% 
#   rename() %>% 
#   column_to_rownames("donorID") %>% 
#   as.matrix()

#### Correlation ####
cor.result <- Hmisc::rcorr(X, type = "pearson")

#### plot ####
corrplot(cor.result$r, p.mat = cor.result$P, method = 'color', type="upper",
         sig.level = c(0.01, 0.05), pch.cex = 0.9, diag=FALSE, 
         insig = 'label_sig', tl.col = "black", col.lim = c(-1,1))

#### Save ####

pdf("publication/FigS5.cell.pct.corr.pdf",
    height = 6, width=6)
corrplot(cor.result$r, p.mat = cor.result$P, method = 'color', type="upper",
         sig.level = c(0.01, 0.05), pch.cex = 0.7, diag=FALSE, 
         insig = 'label_sig', tl.col = "black", tl.cex = 0.7, col.lim = c(-1,1),
         col=colorRampPalette(c("darkblue","white","darkred"))(200))
dev.off()
