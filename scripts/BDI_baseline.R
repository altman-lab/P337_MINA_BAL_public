#For Kelly Colas
library(tidyverse)
library(kimma)
library(limma)
library(RNAetc)

# Data 
load("data_clean/P337_BAL_data.RData")
dat.BAL.abund.norm.voom$targets$age_yrs <- dat.BAL.abund.norm.voom$targets$age_mo/12

# Baseline samples
dat.v4 <- subset_voom(dat.BAL.abund.norm.voom, 
                      lib_filter = "visit == 'V4'")
dat.v5 <- subset_voom(dat.BAL.abund.norm.voom, 
                      lib_filter = "visit == 'V5'")
# Add BDI 
cog <- read_csv("data_raw/addtl.data/summary.scores.final.csv") %>% 
  mutate(visit = case_match(session, 1~"V4", 2~"V5")) %>% 
  mutate(donorID=paste0("MA",subid)) %>% 
  rename(BDI=BDI.total, STAIX=STAIX.total, PSS=PSS.total,
         PANAS_NA=PANAS.NA, PANAS_PA=PANAS.PA)
dat.v4$targets <- dat.v4$targets %>% left_join(cog)
dat.v5$targets <- dat.v5$targets %>% left_join(cog)

# Models
m_BDI4 <- kmFit(dat.v4, patientID="donorID",
                 run_lm = TRUE, model="~BDI",
                 use_weights = TRUE)
m_STAIX4 <- kmFit(dat.v4, patientID="donorID",
               run_lm = TRUE, model="~STAIX",
               use_weights = TRUE)
m_PSS4 <- kmFit(dat.v4, patientID="donorID",
               run_lm = TRUE, model="~PSS",
               use_weights = TRUE)
m_NA4 <- kmFit(dat.v4, patientID="donorID",
               run_lm = TRUE, model="~PANAS_NA",
               use_weights = TRUE)
m_PA4 <- kmFit(dat.v4, patientID="donorID",
              run_lm = TRUE, model="~PANAS_PA",
              use_weights = TRUE)

m_BDI5 <- kmFit(dat.v5, patientID="donorID",
                run_lm = TRUE, model="~BDI",
                use_weights = TRUE)
m_STAIX5 <- kmFit(dat.v5, patientID="donorID",
                  run_lm = TRUE, model="~STAIX",
                  use_weights = TRUE)
m_PSS5 <- kmFit(dat.v5, patientID="donorID",
                run_lm = TRUE, model="~PSS",
                use_weights = TRUE)
m_NA5 <- kmFit(dat.v5, patientID="donorID",
               run_lm = TRUE, model="~PANAS_NA",
               use_weights = TRUE)
m_PA5 <- kmFit(dat.v5, patientID="donorID",
               run_lm = TRUE, model="~PANAS_PA",
               use_weights = TRUE)

save(m_BDI4, m_STAIX4, m_PSS4, m_NA4, m_PA4,
     m_BDI5, m_STAIX5, m_PSS5, m_NA5, m_PA5,
     file="results/other/cognitive_models.RData")

# Results
summarise_kmFit(m_BDI4$lm)
summarise_kmFit(m_NA4$lm)
summarise_kmFit(m_STAIX4$lm)
summarise_kmFit(m_PSS5$lm)
summarise_kmFit(m_NA5$lm)

m_BDI4$lm %>% mutate(visit="Pre-SBPAg") %>% 
  bind_rows(m_NA4$lm %>% mutate(visit="Pre-SBPAg")) %>% 
  bind_rows(m_STAIX4$lm %>% mutate(visit="Pre-SBPAg")) %>% 
  bind_rows(m_PSS5$lm %>% mutate(visit="Post-SBPAg")) %>% 
  bind_rows(m_NA5$lm %>% mutate(visit="Post-SBPAg")) %>% 
  filter(variable != "(Intercept)" & FDR < 0.3) %>% 
  left_join(dat.v4$genes, by=c("gene"="geneName")) %>% 
  select(visit, variable, hgnc_symbol, gene, FDR) %>% 
  write_csv("~/Desktop/P337_MINA_cognitive_DEGs.csv")

# None signif FDR < 0.3
summarise_kmFit(m_PSS4$lm)
summarise_kmFit(m_PA4$lm)
summarise_kmFit(m_BDI5$lm)
summarise_kmFit(m_STAIX5$lm)
summarise_kmFit(m_PA5$lm)
