library(tidyverse)

#### Data ####
##Numeric
meta <- read_csv("data_clean/P337_BAL_meta_clean.csv") %>% 
  distinct(donorID, age_mo, `FEV1 screening WLAC`, `FEV1 PP`) %>% 
  mutate(age_yr = age_mo/12) %>% 
  pivot_longer(-donorID) %>% 
  mutate(visit="none")
  
meta2 <- read_csv("data_clean/P337_BAL_meta_clean.csv") %>% 
  select(donorID, FEV1.pctPP.preAlbuterol_V4:FeNO.PreBro_V5) %>% 
  distinct() %>% 
  pivot_longer(-donorID) %>% 
  separate(name, into=c("name","visit"), sep="_")

meta3 <- read_csv("data_clean/P337_BAL_meta_clean.csv") %>% 
  select(donorID, visit, ends_with(".pct")) %>% 
  distinct() %>% 
  pivot_longer(-c(donorID, visit))

meta4 <- read_csv("data_raw/addtl.data/summary.scores.final.csv")  %>% 
  mutate(donorID=paste0("MA",subid),
         visit = recode(session, `1`="V4",`2`="V5")) %>% 
  select(donorID, visit, BDI.total, ACQ.mean) %>% 
  filter(donorID %in% meta$donorID) %>% 
  pivot_longer(-c(donorID, visit))

meta.num <- bind_rows(meta, meta2, meta3, meta4) %>% 
  arrange(donorID, name)

##Categorical
meta.cat <- read_csv("data_clean/P337_BAL_meta_clean.csv") %>% 
  distinct(donorID, sex, race) %>% 
  pivot_longer(-donorID)

#### Summary stats ###

meta.num %>% 
  group_by(name, visit) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm=TRUE),
            min = min(value, na.rm=TRUE),
            max = max(value, na.rm=TRUE)) %>% 
  arrange(visit) %>% View()

meta.num %>%
  pivot_wider(names_from = visit) %>% 
  select(-none) %>% 
  drop_na() %>% 
  group_by(name) %>% 
  summarise(ttest = t.test(V4, V5, paired = TRUE)["p.value"]$p.value)

meta.cat %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(total = sum(n),
         perc = n/total*100)

