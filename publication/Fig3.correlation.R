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

counts.mod.delta <- counts.mod %>% 
  #separate V4 and V5
  dplyr::select(-libID) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  #shorten module names
  mutate(module=gsub("module_P337_", "", module)) %>% 
  #wide format
  dplyr::select(-V4,-V5) %>% 
  arrange(donorID, module) %>% 
  pivot_wider(names_from = module, values_from = delta)

#### Protein data #### 
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("name","visit"), sep="_") 

plex.delta <- plex %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  drop_na(delta) %>% 
  dplyr::select(-V4,-V5) %>% 
  arrange(name) %>% 
  mutate(name = paste(name, "multiplex", sep=".")) %>% 
  pivot_wider(values_from = delta)

#### IL17A AND EOS MODULE 02 ####
#### plots ####
dat0 <- full_join(counts.mod, plex, by=c("donorID","module"="name", "value", "visit")) %>% 
  full_join(neuro, by=c("donorID","module"="neuro", "value", "visit")) %>% 
  rename(name=module) %>% 
  dplyr::select(-libID) %>% 
  filter(name %in% c("module_P337_BAL_EOS.pct_02", "IL17A", "L_vA_insula")) %>% 
  mutate(name = recode(name, 
                       "module_P337_BAL_EOS.pct_02"="EOS02 module\nlog2 counts per million",
                       "IL17A"="IL17A\nlog10 pg/ml", 
                       "L_vA_insula"="L vA insula"),
         visit = recode_factor(factor(visit), "V4"="Pre","V5"="Post")) %>% 
  group_by(name, donorID) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down")))

plot0 <- dat0 %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="", y="", color="Post - Pre\nchange") +
  facet_wrap(~name, scales="free") +
  scale_color_manual(values=c("down"="#74add1","up"="#d73027"))
plot0

#### corr plots ####
dat <- full_join(counts.mod.delta, plex.delta) %>% 
  full_join(neuro.delta) %>% 
  dplyr::select(donorID, BAL_EOS.pct_02, IL17A.multiplex, L_vA_insula) %>% 
  left_join(dplyr::select(dat0, donorID, diff.col))


corr1 <- rcorr(x=dat$L_vA_insula, y=dat$BAL_EOS.pct_02, type="pearson") 

title1 <- paste(paste("R", signif(corr1$r[1,2], digits=3), sep=" = "),
                paste("P", signif(corr1$P[1,2], digits=3), sep=" = "), sep=", ")

plot1 <- dat %>% 
  drop_na(L_vA_insula, BAL_EOS.pct_02) %>% 
  ggplot(aes(y=L_vA_insula, x=BAL_EOS.pct_02)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE, color="#d73027") +
  labs(y="Post - pre L vA insula",
       x="Post - pre EOS02 module\nlog2 counts per million",
       title=title1) +
  theme_classic()

corr2 <- rcorr(x=dat$L_vA_insula, y=dat$IL17A.multiplex, type="pearson") 

title2 <- paste(paste("R", signif(corr2$r[1,2], digits=3), sep=" = "),
                paste("P", signif(corr2$P[1,2], digits=3), sep=" = "), sep=", ")

plot2 <- dat %>% 
  drop_na(L_vA_insula,IL17A.multiplex) %>% 
  ggplot(aes(y=L_vA_insula, x=IL17A.multiplex)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE, color="#d73027") +
  labs(y="Post - pre L vA insula",
       x="Post - pre IL17A\nlog10 pg/ml",
       title=title2) +
  theme_classic()

#### Save ####
row2 <- plot_grid(plot1,plot2, nrow=1)
row2

ggsave(plot_grid(plot0, row2, nrow=2, labels = c("(A)","(B)")), 
       filename = "publication/Fig3.correlation.png",
       height=6, width=6.5)
ggsave(plot_grid(plot0, row2, nrow=2, labels = c("(A)","(B)")), 
       filename = "publication/Fig3.correlation.pdf",
       height=6, width=6.5)
