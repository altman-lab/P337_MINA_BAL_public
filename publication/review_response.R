library(ppcor)
library(tidyverse)
library(readxl)
library(Hmisc)
library(ggpubr)
library(patchwork)

##### Data #####
#RNAseq
load("data_clean/P337_BAL_data.RData")

#Modules
mod.files <- list.files(path="results/module_level/", pattern="mod_voom_counts",
                        recursive = TRUE, full.names = TRUE)

counts.mod <- data.frame()
for(file in mod.files){
  counts.temp <- read_csv(file) %>% 
    pivot_longer(-c(module,module.char), names_to="libID")
  
  counts.mod <- rbind(counts.mod, counts.temp)
}

#Calculate delta
counts.mod.delta <- counts.mod %>% 
  #Remove modules 00
  filter(!grepl("_00",module)) %>% 
  dplyr::select(-module.char) %>% 
  #Add metadata to denote V4, V5
  full_join(dplyr::select(dat.BAL.abund.norm.voom$targets, 
                          libID, donorID, visit), by = "libID") %>% 
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

#Cytokine protein
plex.delta <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  drop_na(delta) %>% 
  dplyr::select(-V4,-V5) %>% 
  arrange(name) %>% 
  mutate(name = paste(name, "multiplex", sep=".")) %>% 
  pivot_wider(values_from = delta) 

#fMRI
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

##### Combine data #####

dat <- full_join(counts.mod.delta, plex.delta) %>% 
  select(donorID, BAL_EOS.pct_02,
         FLT3L.multiplex,IL17A.multiplex,GCSF.multiplex,IL10.multiplex,PDGFAA.multiplex,
         CCL27.multiplex,IFNA2.multiplex,GMCSF.multiplex,CCL21.multiplex,TGFA.multiplex,
         CCL26.multiplex,IL23.multiplex,CCL7.multiplex,CXCL1.multiplex,IL16.multiplex) %>% 
  full_join(neuro.delta) %>% 
  left_join(dat.BAL.abund.norm.voom$targets %>% 
              mutate(age_yr = age_mo/12) %>% 
              distinct(donorID, age_yr, sex, race) %>% 
              mutate(sex = case_when(sex == "M" ~ 0,
                                     sex == "F" ~ 1),
                     race = case_when(race =="white" ~ 1,
                                      TRUE ~ 1))) %>% 
  arrange(donorID) %>% 
  column_to_rownames("donorID")

#### Correlation of covariates controlled by age/sex ####
# partial correlation between x and y, controlling for z
# pcor.test(x, y, z)
#Run correlations with none, age, or sex covariate

cor.result <- data.frame()

for(i in colnames(dat)[1:24]){
  print(i)
  for(k in colnames(dat)[1:24]){
    temp <- rcorr(dat[,k], dat[,i])
    none <- data.frame(x = k, y = i, covariate = "none", 
                       estimate = temp$r[1,2],
                       p.value = temp$P[1,2])
    dat.temp <- dat[,c(i,k,"age_yr","sex")] %>% drop_na()
    age <- pcor.test(dat.temp[,k], dat.temp[,i], dat.temp$age_yr) %>% 
      mutate(x = k, y = i, covariate = "age", .before = 1)
    sex <- pcor.test(dat.temp[,k], dat.temp[,i], dat.temp$sex) %>% 
      mutate(x = k, y = i, covariate = "sex", .before = 1)
    
    cor.result <- bind_rows(none, age, sex) %>% 
      bind_rows(cor.result)
  }
}

##### Results #####
#filter for comparisons that are significant without any covariates
none_signif <- cor.result %>% 
  filter(covariate=="none" & p.value < 0.05 & x != y) %>% 
  .[!duplicated(data.frame(t(apply(.,1,sort)))),] %>% #remove duplicate that are x,y vs y,x
  distinct(x,y)
dim(none_signif)

#Summarise change in R values
cor.result %>% 
  inner_join(none_signif) %>% 
  select(x,y,covariate,estimate) %>% 
  pivot_wider(names_from = covariate, values_from = estimate) %>% 
  mutate(age_delta = abs(age)-abs(none),
         sex_delta = abs(sex)-abs(none)) %>% 
  select(x,y,contains("delta")) %>% 
  pivot_longer(-c(x:y)) %>% 
  group_by(name) %>% 
  summarise(mean = mean(abs(value)),
            sd = sd(abs(value)),
            abs_min = min(abs(value)),
            abs_max = max(abs(value)))

#How many significant correlations change to non-significant?
pval <- cor.result %>% 
  inner_join(none_signif) %>% 
  select(x,y,covariate,p.value) %>% 
  mutate(p.value = case_when(p.value < 0.05 ~ "Y",
                            TRUE ~ "N")) %>% 
  pivot_wider(names_from = covariate, values_from = p.value) %>% 
  # filter()
  mutate(age_delta = case_when(age==none ~ "no_change", TRUE~age),
         sex_delta = case_when(sex==none ~ "no_change", TRUE~sex)) %>% 
  select(x,y,contains("delta")) %>% 
  pivot_longer(-c(x:y)) 

pval %>% 
  count(name,value)

#Which pairs have changes in significance?
pval %>% 
  filter(value!="no_change") %>% 
  arrange(name) %>% 
  pivot_wider(values_from = y) %>% 
  arrange(x) %>% 
  View()
