# Format and save supplmental tables for publication #
library(tidyverse)
library(readxl)
`%notin%` <- Negate(`%in%`)

#### fMRI ####
#Time point 1 = visit 4 = pre
neuro <- read_excel(sheet="T1",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(visit="pre")

#Time point 2 = visit 5 = post
neuro <- read_excel(sheet="T2",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(visit="post") %>% 
  #Combine with other visit
  full_join(neuro) %>% 
  filter(neuro %notin% c("BDI","LSI","R_insula_3way","LPR_ACC_EOS")) %>% 
  pivot_wider(names_from = neuro) %>% 
  mutate(donorID = paste("MA",idnum,sep="")) %>% 
  select(donorID, visit, everything(), -idnum) %>% 
  arrange(donorID, visit) 

neuro %>% 
  write_csv("publication/TableS2.fMRI.csv")

#### Log2 CPM genes ####
attach("data_clean/P337_BAL_data.RData")

as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  #Add HGNC names
  rownames_to_column("geneName") %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  select(geneName,hgnc_symbol, everything(), -gene_biotype) %>% 
  rename(ENSEMBL=geneName) %>% 
  #Save
  write_csv("publication/TableS4.gene.counts.csv")

#### Log2 CPM modules ####
attach("data_clean/P337_BAL_module_data.RData")

mod.voom %>% 
  #Clean module names
  mutate(module = gsub("P337_","",module),
         module = gsub(".pct","",module)) %>% 
  #Save
  write_csv("publication/TableS7.module.counts.csv")

#### Library metadata ####
dat.BAL.abund.norm.voom$targets %>% 
  select(libID, norm.factors, 
         donorID, visit, group, flow_cell,
         total_sequences, median_cv_coverage, mapped_reads_w_dups, 
         EOS.pct, PMN.pct) %>% 
  rename(type=group, batch=flow_cell, TMM_norm_factor=norm.factors,
         EOS_percent=EOS.pct, PMN_percent=PMN.pct) %>% 
  mutate(visit = recode(visit, "V4"="pre", "V5"="post")) %>% 
  write_csv("publication/TableS3.RNAseq.library.metadata.csv")

#### Gene linear models ####
gene.visit <- read_csv("results/gene_level/P337_BAL_gene_visit.csv") %>% 
  filter(group == "visit" & adj.P.Val <= 0.3 ) %>% 
  #Rename specific to this model
  select(geneName, logFC, adj.P.Val) %>% 
  rename(visit_logFC=logFC, visit_FDR=adj.P.Val)

gene.EOS <- read_csv("results/gene_level/P337_BAL_gene_EOS.csv") %>% 
  filter(group == "EOS.pct") %>% 
  #Rename specific to this model
  select(geneName, logFC, adj.P.Val) %>% 
  rename(EOS_logFC=logFC, EOS_FDR=adj.P.Val)

gene.PMN <- read_csv("results/gene_level/P337_BAL_gene_PMN.csv") %>% 
  filter(group == "PMN.pct") %>% 
  #Rename specific to this model
  select(geneName, logFC, adj.P.Val) %>% 
  rename(PMN_logFC=logFC, PMN_FDR=adj.P.Val)

#Module assignment
mod <- mod.genes %>% 
  select(geneName, hgnc_symbol, module) %>% 
  #Clean module names
  mutate(module = gsub("P337_","",module),
         module = gsub(".pct","",module)) %>% 
  #separate EOS and PMN modules
  mutate(name = ifelse(grepl("EOS", module), "EOS_module",
                        ifelse(grepl("PMN", module), "PMN_module", NA))) %>% 
  pivot_wider(values_from = module)

full_join(gene.visit, gene.EOS) %>% 
  full_join(gene.PMN) %>% 
  #add module assignment and HGNC
  full_join(mod) %>% 
  rename(ENSEMBL=geneName) %>% 
  select(ENSEMBL, hgnc_symbol, everything()) %>% 
  write_csv("publication/TableS6.gene.linear.models.csv")

#### Correlation Post-Pre ####
#Add if selected by SPLS
load("results/PLS/SPLS.RData")
splsX <- as.data.frame(result.spls$loadings$X) %>% 
  filter(comp1 != 0 | comp2 != 0)

read_csv("results/PLS/pearson.correlation.csv") %>% 
  mutate(SPLS_selected = ifelse(X_variable %in% rownames(splsX) |
                                  X_variable == "module_BAL_EOS_02", "Y",NA)) %>% 
  dplyr::select(X_variable, SPLS_selected, everything()) %>% 
  arrange(SPLS_selected) %>% 
  write_csv("publication/TableS8.Pearson.correlation.csv")

#### ELISA ####

read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("name","visit")) %>% 
  mutate(visit = recode(visit, "V4"="pre","V5"="post")) %>% 
  rename(donorID=ptID) %>% 
  pivot_wider() %>% 
  arrange(donorID, desc(visit)) %>% 
  write_csv("publication/TableS5.ELISA.csv")
