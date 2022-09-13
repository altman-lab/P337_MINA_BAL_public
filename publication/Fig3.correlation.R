library(tidyverse)
library(readxl)
library(cowplot)
library(Hmisc)

#### Neuro data ####
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

#### Gene module data ####
attach("data_clean/P337_BAL_data.RData")

counts.mod <- read_csv("results/module_level//module_P337_BAL_EOS.pct_deepSplit2_minMod50/P337_BAL_EOS.pct_mod_voom_counts.csv") %>% 
  pivot_longer(-c(module,module.char), names_to="libID") %>% 
  #Remove modules 00
  filter(!grepl("_00",module)) %>% 
  dplyr::select(-module.char) %>% 
  #Add metadata to denote V4, V5
  full_join(dplyr::select(dat.BAL.abund.norm.voom$targets, 
                          libID, donorID, visit), by = "libID")

#### Protein data #### 
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("name","visit"), sep="_") 

#### Cell data ####
load("data_clean/P337_BAL_data.RData")

cell <- dat.BAL.abund.norm.voom$targets %>% 
  select(visit, donorID, EOS.pct, NEUT.pct) %>% 
  pivot_longer(EOS.pct:NEUT.pct)

#### IL17A AND EOS MODULE 02 ####
#### plots ####
OI <-  c("module_P337_BAL_EOS.pct_02", "IL17A", "EOS.pct", "L_vA_insula")

dat0 <- full_join(counts.mod, plex, by=c("donorID","module"="name", "value", "visit")) %>% 
  full_join(neuro, by=c("donorID","module"="neuro", "value", "visit")) %>% 
  rename(name=module) %>% 
  full_join(rownames_to_column(cell, "libID")) %>% 
  dplyr::select(-libID) %>% 
  filter(name %in% OI) %>% 
  #remove failed fMRI
  filter(donorID != "MA1012") %>% 
  mutate(name = recode_factor(factor(name), 
                              "module_P337_BAL_EOS.pct_02"="EOS02 module\nlog2 CPM",
                              "IL17A"="IL-17A\nlog10 pg/ml", 
                              "EOS.pct"="Eosinophil\n cell %",
                              "L_vA_insula"="L vA insula\nBOLD score"),
         visit = recode_factor(factor(visit), "V4"="Pre","V5"="Post")) %>% 
  group_by(name, donorID) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) 

plot00 <- dat0 %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmental bronchial provocation with allergen (SBP-Ag)",
       y="", color="Post - Pre\nchange") +
  facet_wrap(~name, scales="free", nrow=1) +
  scale_color_manual(values=c("down"="#74add1","up"="#d73027"))
# plot00

#### corr plots ####
dat.delta <- dat0 %>% 
  distinct(name, donorID, diff) %>% 
  pivot_wider(values_from = diff)
  
corr0 <- rcorr(x=dat.delta$`L vA insula\nBOLD score`,
               y=dat.delta$`Eosinophil\n cell %`, type="pearson") 

title0 <- paste(paste("R", signif(corr0$r[1,2], digits=2), sep=" = "),
                paste("P", signif(corr0$P[1,2], digits=2), sep=" = "), sep=", ")

plot0 <- dat.delta %>% 
  ggplot(aes(y=`L vA insula\nBOLD score`, x=`Eosinophil\n cell %`)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE, color="#d73027") +
  labs(y="Post - pre L vA insula",
       x="Post - pre eosinophil\ncell %",
       title=title0) +
  theme_classic()

corr1 <- rcorr(x=dat.delta$`L vA insula\nBOLD score`, 
               y=dat.delta$`EOS02 module\nlog2 CPM`, type="pearson") 

title1 <- paste(paste("R", signif(corr1$r[1,2], digits=2), sep=" = "),
                paste("P", signif(corr1$P[1,2], digits=2), sep=" = "), sep=", ")

plot1 <- dat.delta %>% 
  ggplot(aes(y=`L vA insula\nBOLD score`, x=`EOS02 module\nlog2 CPM`)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE, color="#d73027") +
  labs(y="Post - pre L vA insula",
       x="Post - pre EOS02 module\nlog2 counts per million",
       title=title1) +
  theme_classic()

corr2 <- rcorr(x=dat.delta$`L vA insula\nBOLD score`,
               y=dat.delta$`IL-17A\nlog10 pg/ml`, type="pearson") 

title2 <- paste(paste("R", signif(corr2$r[1,2], digits=2), sep=" = "),
                paste("P", signif(corr2$P[1,2], digits=2), sep=" = "), sep=", ")

plot2 <- dat.delta %>% 
  ggplot(aes(y=`L vA insula\nBOLD score`, x=`IL-17A\nlog10 pg/ml`)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE, color="#d73027") +
  labs(y="Post - pre L vA insula",
       x="Post - pre IL-17A\nlog10 pg/ml",
       title=title2) +
  theme_classic()

#### Save ####
row2 <- plot_grid(plot1,plot2,plot0, nrow=1)
row2

ggsave(plot_grid(plot00, NULL, row2, nrow=3, labels = c("(A)","","(B)"),
                 rel_heights = c(1, 0.05, 1)), 
       filename = "publication/Fig3.correlation.png",
       height=5, width=7)
ggsave(plot_grid(plot00, NULL, row2, nrow=3, labels = c("(A)","","(B)"),
                 rel_heights = c(1, 0.05, 1)), 
       filename = "publication/Fig3.correlation.pdf",
       height=5, width=7)
