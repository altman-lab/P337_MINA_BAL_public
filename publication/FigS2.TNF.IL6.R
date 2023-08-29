library(tidyverse)
library(readxl)
library(patchwork)
library(ggpubr)
library(Hmisc)

#### Targeted gene-protein pairs ####
genes <- c("TPSAB1", "CPA3", "KIT","IL6", "IL13")
prots <- c("TNFA", "IL6", "RANTES", "IL8", "EGF", "IL5","IL13")

genes.OI <- c(paste(genes, "gene"), paste(prots, "protein"))

#### Data ####
# Gene data
load("data_clean/P337_BAL_data.RData")

# Module data
load("data_clean/P337_BAL_module_data.RData")

# Protein data
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  rename(hgnc_symbol=name)

# fMRI data
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

#### Stats ####
gene_lm <- read_csv("results/gene_level/P337_BAL_gene_visit.csv") %>% 
  filter(group == "visit") %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  filter(hgnc_symbol %in% genes) %>% 
  select(hgnc_symbol, adj.P.Val) %>% 
  mutate(group="gene") %>% 
  mutate(adj.P.Val = formatC(adj.P.Val, digits=1, format="e"))

mod_lm <- read_csv("results/module_level/P337_BAL_mod_visit.csv") %>% 
  filter(group == "visit") %>% 
  filter(geneName %in% c("module_P337_BAL_EOS.pct_01")) %>% 
  select(geneName, adj.P.Val) %>% 
  mutate(group="module") %>% 
  mutate(adj.P.Val = formatC(adj.P.Val, digits=1, format="e"))

prot_lm <- read_csv("results/protein_level/P337_BAL_protein_visit.csv")%>% 
  filter(group == "visit") %>% 
  filter(geneName %in% prots) %>% 
  rename(hgnc_symbol = geneName) %>% 
  select(hgnc_symbol, adj.P.Val)%>% 
  mutate(group="protein") %>% 
  mutate(adj.P.Val = formatC(adj.P.Val, digits=1, format="e"))

#### Combine ####
dat_gene <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  #gene expression
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "libID") %>% 
  #sample metadata
  left_join(dat.BAL.abund.norm.voom$targets %>% 
              select(libID,donorID,visit)) %>% 
  #gene names
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  filter(hgnc_symbol %in% genes) %>% 
  #add gene stats
  left_join(gene_lm) %>% 
  #fix labels
  mutate(visit=recode_factor(visit, "V4"="Pre", "V5"="Post"))  %>% 
  group_by(hgnc_symbol, donorID) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) 

dat_mod <- mod.voom %>% 
  #gene expression
  pivot_longer(-module, names_to = "libID") %>% 
  #sample metadata
  left_join(dat.BAL.abund.norm.voom$targets %>% 
              select(libID,donorID,visit)) %>% 
  filter(module %in% c("module_P337_BAL_EOS.pct_01")) %>% 
  #add gene stats
  left_join(mod_lm, by=c("module"="geneName")) %>% 
  #fix labels
  mutate(visit=recode_factor(visit, "V4"="Pre", "V5"="Post"))  %>% 
  group_by(module, donorID) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) %>% 
  mutate(module = gsub("module_P337_BAL_|.pct|_", "", module))

dat_prot <- plex %>% 
  filter(hgnc_symbol %in% prots) %>% 
  #add prot stats
  left_join(prot_lm) %>% 
  #fix labels
  mutate(visit=recode_factor(visit, "V4"="Pre", "V5"="Post")) %>% 
  group_by(hgnc_symbol, donorID) %>% 
  mutate(diff = value[visit=="Post"]-value[visit=="Pre"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) 

#### Pre vs post plots ####
p_prot1 <- dat_prot %>% 
  mutate(facet.lab = paste0(hgnc_symbol, " protein\nFDR = ", adj.P.Val)) %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmental bronchial provocation with allergen (SBP-Ag)",
       y="Protein log10 pg/ml", color="Post - Pre\nchange") +
  facet_wrap(~facet.lab, scales="free", nrow=1) +
  scale_color_manual(values=c("down"="#74add1","up"="#FF8679")) 
# p_prot1

p_gene1 <- dat_gene %>% 
  mutate(facet.lab = paste0(hgnc_symbol, " gene\nFDR = ", adj.P.Val)) %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmental bronchial provocation with allergen (SBP-Ag)",
       y="Gene log2 CPM", color="Post - Pre\nchange") +
  facet_wrap(~facet.lab, scales="free", nrow=1) +
  scale_color_manual(values=c("down"="#74add1","up"="#FF8679"))
# p_gene1

p_mod1 <- dat_mod %>% 
  mutate(facet.lab = paste0(module, " module\nFDR = ", adj.P.Val)) %>% 
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmental bronchial\nprovocation with allergen\n(SBP-Ag)",
       y="Module mean log2 CPM", color="Post - Pre\nchange") +
  facet_wrap(~facet.lab, scales="free", nrow=1) +
  scale_color_manual(values=c("down"="#74add1","up"="#FF8679")) +
  theme(legend.position = "none")
# p_mod1

#### delta corr to L vA insula ####
#proteins
plot_ls <- list()

for (p in sort(unique(dat_prot$hgnc_symbol))){
  prot_temp <- dat_prot %>% 
    filter(hgnc_symbol==p) %>% 
    pivot_wider(names_from = visit) %>% 
    mutate(delta=Post-Pre)  %>% 
    select(-Pre,-Post) %>% 
    left_join(neuro.delta %>% select(donorID, L_vA_insula)) %>% 
    drop_na(L_vA_insula, delta)
  
  corr0 <- rcorr(x=prot_temp$L_vA_insula,
                 y=prot_temp$delta, type="pearson") 
  
  title0 <- paste(paste("R", signif(corr0$r[1,2], digits=2), sep=" = "),
                  paste("P", signif(corr0$P[1,2], digits=2), sep=" = "), sep=", ")
  
  if(p==sort(unique(dat_prot$hgnc_symbol))[1]){
    ylab <-"Post - pre L vA insula"
  } else{
    ylab <- ""
  }
  
  plot_ls[[p]] <- prot_temp %>% 
    ggplot(aes(y=L_vA_insula, x=delta)) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE, color="#FF8679", formula = 'y ~ x') +
    labs(y=ylab,
         x=paste("Post - pre", p, "\nlog10 pg/ml"),
         title=title0) +
    theme_classic()+
    theme(plot.title  = element_text(size=10))
}

#genes
plot_ls2 <- list()

for (g in sort(unique(dat_gene$hgnc_symbol))){
  gene_temp <- dat_gene %>% 
    filter(hgnc_symbol==g) %>% 
    select(-libID) %>% 
    pivot_wider(names_from = visit) %>% 
    mutate(delta=Post-Pre)  %>% 
    select(-Pre,-Post) %>% 
    left_join(neuro.delta %>% select(donorID, L_vA_insula)) %>% 
    drop_na(L_vA_insula, delta)
  
  corr0 <- rcorr(x=gene_temp$L_vA_insula,
                 y=gene_temp$delta, type="pearson") 
  
  title0 <- paste(paste("R", signif(corr0$r[1,2], digits=2), sep=" = "),
                  paste("P", signif(corr0$P[1,2], digits=2), sep=" = "), sep=", ")
  
  if(g==sort(unique(dat_gene$hgnc_symbol))[1]){
    ylab <-"Post - pre L vA insula"
  } else{
    ylab <- ""
  }
  
  plot_ls2[[g]] <- gene_temp %>% 
    ggplot(aes(y=L_vA_insula, x=delta)) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE, color="#FF8679", formula = 'y ~ x') +
    labs(y=ylab,
         x=paste("Post - pre", g, "\nlog2 CPM"),
         title=title0) +
    theme_classic()+
    theme(plot.title  = element_text(size=10))
}

#modules
plot_ls3 <- list()

for (m in sort(unique(dat_mod$module))){
  mod_temp <- dat_mod %>% 
    filter(module==m) %>% 
    select(-libID) %>% 
    pivot_wider(names_from = visit) %>% 
    mutate(delta=Post-Pre)  %>% 
    select(-Pre,-Post) %>% 
    left_join(neuro.delta %>% select(donorID, L_vA_insula)) %>% 
    drop_na(L_vA_insula, delta)
  
  corr0 <- rcorr(x=gene_temp$L_vA_insula,
                 y=gene_temp$delta, type="pearson") 
  
  title0 <- paste(paste("R", signif(corr0$r[1,2], digits=2), sep=" = "),
                  paste("P", signif(corr0$P[1,2], digits=2), sep=" = "), sep=", ")
  
  if(m==sort(unique(dat_mod$module))[1]){
    ylab <-"Post - pre L vA insula"
  } else{
    ylab <- ""
  }
  
  plot_ls3[[m]] <- mod_temp %>% 
    ggplot(aes(y=L_vA_insula, x=delta)) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE, color="#FF8679", formula = 'y ~ x') +
    labs(y=ylab,
         x=paste("Post - pre", m, "\nlog2 CPM"),
         title=title0) +
    theme_classic()+
    theme(plot.title  = element_text(size=10))
}

corr_plots1 <- wrap_plots(plot_ls, nrow = 1, tag_level="new")
corr_plots2 <- wrap_plots(plot_ls2, nrow = 1, tag_level="new")
corr_plots3 <- wrap_plots(plot_ls3, nrow = 1, tag_level="new")


#### Save ####
lo <- "
ACCCCC
DFFFFF
BBBBBB
EEEEEE
"
plot_all <- 
  p_mod1 + p_prot1 + p_gene1 + corr_plots3 + corr_plots1 + corr_plots2 + 
  plot_layout(design = lo) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

# plot_all

ggsave(plot_all, filename="publication/FigS2.TNF.IL6.png",
       width=14, height=10)
ggsave(plot_all, filename="publication/FigS2.TNF.IL6.pdf",
       width=14, height=10)
