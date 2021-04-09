library(tidyverse)

#### Data ####
meta <- read_csv("data_clean/P337_BAL_meta_clean.csv") %>% 
  distinct(donorID, age_mo, sex, `FEV1 screening WLAC`, `FEV1 PP`)

#### Summary stats ###

meta %>% 
  summarise(across(c(age_mo, `FEV1 screening WLAC`, `FEV1 PP`), ~mean(.)))
meta %>% 
  summarise(across(c(age_mo, `FEV1 screening WLAC`, `FEV1 PP`), ~sd(.)))

table(meta$sex)
