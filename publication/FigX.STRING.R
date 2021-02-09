library(tidyverse)
source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/STRING_network_fxn.R")
set.seed(3698)

#### Data modules ####
#Genes in modules.
load("data_clean/P337_BAL_module_data.RData")

#Format Hallmark terms for coloring.
H.format <- read_csv("results/enrichment/enrich_modules_H.csv") %>% 
  #significant
  filter(size.overlap.term >= 2 & p.adjust <= 0.3) %>% 
  #Beautify labels
  mutate(Description = gsub("HALLMARK_","",Description),
         Description = gsub("_"," ", Description))

#### Plot genes in modules ####
for(mod in c("EOS.pct_02")){
  mod.genes.temp <- mod.genes %>% 
    filter(grepl(mod, module)) %>% 
    distinct(hgnc_symbol) %>% unlist(use.names=FALSE)
  
  H.format.temp <- H.format %>% 
    filter(grepl(mod, group))
  
  string.plot(genes=mod.genes.temp, version="11", score_threshold=700,
              layout='fr', discard='edge.keep.enrich',
              enrichment=H.format.temp, ID="SYMBOLs",
              outdir="publication/", basename=paste("publication/FigX2.", mod, ".", sep=""),
              width=8, height=8)
}

#### Data SPLS ####
load("results/PLS/SPLS.RData")

cyto <- as.data.frame(result.spls$loadings$X) %>% 
  rownames_to_column("var") %>% 
  filter(comp1 != 0 | comp2 != 0) %>% 
  filter(!grepl("BAL EOS", var)) %>% 
  mutate(gene = recode(gsub(".multiplex","",var),
                       "FLT3L"="FLT3LG",
                       "GCSF"="CSF3",
                       "GMCSF"="CSF2",
                       "IL12p70"="IL12A_IL12B",
                       "IL1RA"="IL1RN",
                       "IL23"="IL23A",
                       "IL8"="CXCL8",
                       "IP10"="CXCL10",
                       "PDGFAA"="PDGFA",
                       "PDGFAB"="PDGFA_PDGFB",
                       "RANTES"="CCL5",
                       "sCD40L"="CD40LG",
                       "TARC"="CCL17",
                       "TNFA"="TNF",
                       "TNFB"="LTA",
                       "VEGF"="VEGFA_VEGFB_VEGFC")) %>% 
  #Unnest multi-annotations
  separate(gene, into=c("a","b","c"), sep="_") %>% 
  pivot_longer(a:c, values_to = "gene") %>% 
  drop_na(gene)

#### Plot SPLS cytokines ####
string.plot(genes=cyto$gene, version="11", score_threshold=700,
            layout='fr', 
            ID="SYMBOLs", 
            outdir="publication/", basename="FigX.cyto700",
            width=8, height=8)
