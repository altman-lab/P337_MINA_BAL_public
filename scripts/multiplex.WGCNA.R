library(tidyverse)
set.seed(4389)

#### Data ####
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-ptID) %>% 
  #Log10 transform
  mutate(value = ifelse(value==0,0, log10(value))) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  pivot_wider() %>% 
  mutate(rowname=paste(ptID,visit, sep="_")) %>% 
  select(-ptID,-visit) %>% 
  column_to_rownames()

#### Build modules ####
library(WGCNA)

#Calculate thresholds
sft <- pickSoftThreshold(plex, 
                         powerVector=c(1:10), verbose=5,
                         networkType = "signed")
#Select threshold 
power.t <- 2

    
mod.net <- blockwiseModules(plex,
                            power=power.t, 
                            networkType="signed",
                            TOMType="signed",
                            corType="bicor",
                            maxBlockSize=100,
                            minModuleSize=5,
                            deepSplit=4, 
                            numericLabels=TRUE,
                            saveTOMFileBase="TOM-blockwise",
                            nthreads=4, impute=FALSE)
#Extract results  
mods <- as.data.frame(mod.net$colors) %>% 
  rownames_to_column("protein") %>% 
  dplyr::rename(module = "mod.net$colors") %>% 
  #add leading 0 to module names for correct sorting of factor
  mutate(module.char = ifelse(module <= 9, 
                              paste("0", module, sep=""),
                              module)) %>% 
  #Add color var
  mutate(mod.color = labels2colors(mod.net$colors))

#Number of genes in each module
table(mods$module.char)

#### Delta ####
plex.delta <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-ptID) %>% 
  #Log10 transform
  mutate(value = ifelse(value==0,0, log10(value))) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  drop_na(delta) %>% 
  select(-V4,-V5) %>% 
  arrange(name) %>% 
  mutate(name = paste(name, "multiplex", sep=".")) %>% 
  pivot_wider(values_from = delta) %>% 
  column_to_rownames("ptID")

#Calculate thresholds
sft <- pickSoftThreshold(plex.delta, 
                         powerVector=c(1:10), verbose=5,
                         networkType = "signed")
#Select threshold 
power.t <- 2


mod.net <- blockwiseModules(plex.delta,
                            power=power.t, 
                            networkType="signed",
                            TOMType="signed",
                            corType="pearson",
                            maxBlockSize=100,
                            minModuleSize=5,
                            deepSplit=4, 
                            numericLabels=TRUE,
                            saveTOMFileBase="TOM-blockwise",
                            nthreads=4, impute=FALSE)
#Extract results  
mods <- as.data.frame(mod.net$colors) %>% 
  rownames_to_column("protein") %>% 
  dplyr::rename(module = "mod.net$colors") %>% 
  #add leading 0 to module names for correct sorting of factor
  mutate(module.char = ifelse(module <= 9, 
                              paste("0", module, sep=""),
                              module)) %>% 
  #Add color var
  mutate(mod.color = labels2colors(mod.net$colors))

#Number of genes in each module
table(mods$module.char)

#### Correlation with RNAseq ####
#### Delta RNAseq modules ####
#Read in all module expression
mod.files <- list.files(path="results/module_level/", pattern="mod_voom_counts",
                        recursive = TRUE, full.names = TRUE)

counts.all <- data.frame()
for(file in mod.files){
  counts.temp <- read_csv(file) %>% 
    pivot_longer(-module, names_to="libID")
  
  counts.all <- rbind(counts.all, counts.temp)
}

#Metadata
load("data_clean/P337_BAL_data.RData")
load("data_clean/P337_BE_data.RData")

meta.all <- bind_rows(dat.BAL.abund.norm.voom$targets, 
                      dat.BE.abund.norm.voom$targets)

#Calculate delta
counts.delta <- counts.all %>% 
  #Remove modules 00
  filter(!grepl("_00",module)) %>% 
  #Add metadata to denote V4, V5
  full_join(select(meta.all, libID, donorID, visit), by = "libID") %>% 
  #separate V4 and V5
  select(-libID) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  #shorten module names
  mutate(module=gsub("module_P337_", "", module)) %>% 
  #wide format
  select(-V4,-V5) %>% 
  arrange(donorID, module) %>% 
  pivot_wider(names_from = module, values_from = delta)

#### Delta multiplex modules ####
## Means
plex.delta.mod <- plex.delta %>% 
  rownames_to_column("sampID") %>% 
  pivot_longer(-sampID, names_to = "protein") %>% 
  inner_join(mods) %>% 
  filter(module != 0) %>% 
  group_by(sampID, module.char) %>% 
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>% 
  mutate(module.char=paste("module_multiplex", module.char, sep="_")) %>% 
  pivot_wider(names_from = module.char, values_from = mean)

write_csv(plex.delta.mod, "plex.delta.mod.csv")
#### Corr ####
delta.all <- full_join(counts.delta, plex.delta.mod, by=c("donorID"="sampID")) %>% 
  select(!contains("BE_")) %>% 
  column_to_rownames("donorID")

library(Hmisc)
library(corrplot)

corr.result <- rcorr(as.matrix(delta.all), type="pearson")
  corr.R <- corr.result$r
  corr.P <- corr.result$P

#### Calculate FDR 
P.lower <- corr.result$P
  P.lower[upper.tri(P.lower, diag=FALSE)] <- NA
  
  FDR.lower <-  P.lower %>%
    rownames_to_column("variable") %>% 
    pivot_longer(-variable) %>% 
    #Variables of interest
    filter(variable %in% all_vars & name %in% all_vars) %>% 
    #FDR groups
    mutate(value = p.adjust(value, method="BH")) 
  
  P.upper <-  P %>% 
    column_to_rownames("variable")
  P.upper[lower.tri(P.upper, diag=FALSE)] <- NA
  
  FDR.upper <-  P.upper %>%
    rownames_to_column("variable") %>% 
    pivot_longer(-variable) %>% 
    filter(variable %in% all_vars & name %in% all_vars) %>% 
    mutate(value = p.adjust(value, method="BH"))
  
  FDR <- bind_rows(FDR.lower, FDR.upper) %>% 
    distinct() %>% 
    drop_na(value) %>% 
    pivot_wider()
  
corrplot(as.matrix(corr.R), method="color", order="original",
           #Reverse color order from default
           col=colorRampPalette(c("darkblue", "white",
                                  "darkred"))(20),
           #Change labels
           tl.col="black", tl.srt=90, cl.lim=c(-1,1),
           #Change colorlegend
           cl.pos="b", cl.length = 5,
           #Add significant labels
           p.mat = as.matrix(corr.P), insig = "label_sig",
           sig.level = c(0.01, 0.05), 
           pch.cex=0.9, pch.col="white",
           #Add title
           title="P *<0.05 **<0.01", 
           mar = c(0,0,2,0))
  