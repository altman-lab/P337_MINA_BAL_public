# Format and save supplmental tables for publication #
library(tidyverse)
library(readxl)
library(openxlsx)
`%notin%` <- Negate(`%in%`)

#### Patient metadata ####
meta1 <- read_excel("data_raw/addtl.data/2022-08_Patient_metadata.xlsx") %>% 
  rename(donorID=id, age_months=age.months, sex=Sex, bmi=BMI, 
         FEV1_PP=`Pre-albuterol_FEV1_%_predicted`,
         `FEV1/FVC`=`Pre-albuterol_FEV1/FVC`) %>% 
  select(-age.years)

#### Library metadata ####
meta2 <- dat.BAL.abund.norm.voom$targets %>% 
  select(libID, norm.factors, 
         donorID, visit, group, flow_cell,
         total_sequences, median_cv_coverage, mapped_reads_w_dups, 
         EOS.pct, NEUT.pct) %>% 
  rename(type=group, batch=flow_cell, TMM_norm_factor=norm.factors,
         EOS_percent=EOS.pct, NEUT_percent=NEUT.pct) %>% 
  mutate(visit = recode(visit, "V4"="pre", "V5"="post"))

#### Table S1 ####
write.xlsx(list("patient"=meta1, "rnaseq"=meta2), 
           file = "publication/TableS1.metadata.xlsx")

#### Log2 CPM genes ####
attach("data_clean/P337_BAL_data.RData")

gene <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  #Add HGNC names
  rownames_to_column("geneName") %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  select(geneName,hgnc_symbol, everything(), -gene_biotype) %>% 
  rename(ENSEMBL=geneName)

#### Log10 ELISA ####
elisa <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  pivot_longer(-ptID) %>% 
  mutate(value = log10(value)) %>% 
  separate(name, into=c("name","visit")) %>% 
  mutate(visit = recode(visit, "V4"="pre","V5"="post")) %>% 
  rename(donorID=ptID) %>% 
  pivot_wider() %>% 
  arrange(donorID, desc(visit)) 

#### Log2 CPM modules ####
attach("data_clean/P337_BAL_module_data.RData")

mod <- mod.voom %>% 
  #Clean module names
  mutate(module = gsub("P337_","",module),
         module = gsub(".pct","",module))

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

##### Table S2 ####
write.xlsx(list("gene"=gene, "protein"=elisa, "module"=mod, "fMRI"=neuro), 
           file = "publication/TableS2.data.xlsx")

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

gene.NEUT <- read_csv("results/gene_level/P337_BAL_gene_NEUT.csv") %>% 
  filter(group == "NEUT.pct") %>% 
  #Rename specific to this model
  select(geneName, logFC, adj.P.Val) %>% 
  rename(NEUT_logFC=logFC, NEUT_FDR=adj.P.Val)

#Module assignment
mod.key <- mod.genes %>% 
  select(geneName, hgnc_symbol, module) %>% 
  #Clean module names
  mutate(module = gsub("P337_","",module),
         module = gsub(".pct","",module)) %>% 
  #separate EOS and NEUT modules
  mutate(name = ifelse(grepl("EOS", module), "EOS_module",
                       ifelse(grepl("NEUT", module), "NEUT_module", NA))) %>% 
  pivot_wider(values_from = module)

model.result <- full_join(gene.visit, gene.EOS) %>% 
  full_join(gene.NEUT) %>% 
  #add module assignment and HGNC
  full_join(mod.key) %>% 
  rename(ENSEMBL=geneName) %>% 
  select(ENSEMBL, hgnc_symbol, everything())

#simplfy genes in modules lists
temp <- model.result %>% 
  select(hgnc_symbol, EOS_module, NEUT_module) %>% 
  pivot_longer(-hgnc_symbol) %>% 
  drop_na(value)

mod.ls <- list()
for(mod in sort(unique(temp$value))){
  mod.ls[[mod]] <- temp %>% 
    filter(value == mod) %>% 
    arrange(hgnc_symbol) %>% 
    pull(hgnc_symbol)
}

gene.in.mod <- plyr::ldply(mod.ls, rbind) %>% 
  column_to_rownames(".id") %>% 
  t() %>% as.data.frame()

#### Enrichment ####
enrich <- filter(read_csv("results/enrichment/enrich_sPLS_H.csv"), 
       FDR<=0.2) %>% 
  bind_rows(filter(read_csv("results/enrichment/enrich_sPLS_C2_CP.csv"),
                   FDR<=0.05)) %>% 
  bind_rows(filter(read_csv("results/enrichment/enrich_sPLS_C5_GO.BP.csv"), 
                   FDR<=0.05)) %>% 
  bind_rows(filter(read_csv("results/enrichment/enrich_sPLS_C7.csv"), 
                   FDR<=0.03)) %>% 
  arrange(gs_cat, gs_subcat, FDR) %>% 
  dplyr::select(gs_cat, gs_subcat, pathway, group_in_pathway, 
         size_pathway, pval, FDR, genes) %>% 
  rename(`SPLS genes in term (k)`=group_in_pathway, `total genes in term (K)`=size_pathway,
         category=gs_cat, subcategory=gs_subcat, term=pathway)

#### Correlation Post-Pre ####
#Add if selected by SPLS
load("results/PLS/SPLS.RData")
splsX <- as.data.frame(result.spls$loadings$X) %>% 
  filter(comp1 != 0 | comp2 != 0)

corr <- read_csv("results/PLS/pearson.correlation.csv") %>% 
  mutate(SPLS_selected = ifelse(X_variable %in% rownames(splsX) |
                                  X_variable == "module_BAL_EOS_02", "Y",NA)) %>% 
  dplyr::select(X_variable, SPLS_selected, everything()) %>% 
  arrange(SPLS_selected)

#### Table S3 ####
write.xlsx(list("linear_model" = model.result, "genes_in_modules" = gene.in.mod,
                "enrichment"=enrich, "correlation"=corr),
           file = "publication/TableS3.stats.xlsx")
