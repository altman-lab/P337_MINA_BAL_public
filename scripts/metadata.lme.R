library(tidyverse)
library(broom)
library(lme4)
library(car)

#### Cell % ####
cell.pct <- read_csv("data_clean/P337_BAL_cell_pcts.csv") %>% 
  rename(donorID=ptID)

results <- data.frame()
for(cell in colnames(cell.pct)[3:7]){
  temp <- tidy(Anova(lmer(get(cell) ~ visit + (1|donorID), data=cell.pct))) %>% 
    mutate(cell=cell)
  
  results <- bind_rows(results, temp)
}

results <- results %>% 
  mutate(FDR = p.adjust(p.value))

fc <- cell.pct %>% 
  group_by(visit) %>% 
  summarise(across(EOS.pct:Epi.pct, ~mean(.,na.rm=TRUE))) %>% 
  ungroup() %>% 
  column_to_rownames("visit") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(V5-V4)

#### Multiplex ####  
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  pivot_wider() %>% 
  #Log10 transform
  mutate(across(IL5:CCL15, ~ifelse(.==0,0, log10(.))))

results <- data.frame()
for(prot in colnames(plex)[3:58]){
  temp <- tidy(Anova(lmer(get(prot) ~ visit + (1|donorID), data=plex))) %>% 
    mutate(prot=prot)
  
  results <- bind_rows(results, temp)
}

results <- results %>% 
  mutate(FDR = p.adjust(p.value))

fc <- plex %>% 
  group_by(visit) %>% 
  summarise(across(IL5:CCL15, ~mean(.,na.rm=TRUE))) %>% 
  ungroup() %>% 
  column_to_rownames("visit") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(V5-V4) %>% 
  mutate(unlog = 10^`V5 - V4`)
