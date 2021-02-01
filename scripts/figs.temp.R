library(tidyverse)
library(cowplot)
`%notin%` <- Negate(`%in%`)

#### Cell percentages ####
cell.count <- read_csv("data_raw/addtl.data/P337_cell.counts.csv")

plot1 <- cell.count %>% 
  #convert BAL to number per ml
  mutate(across(c(BAL_EOS_no_V4, BAL_LYM_no_V4, BAL_PMN_no_V4, 
                  BAL_MONO_no_V4), ~./BAL_volume.ml_V4/1E3)) %>% 
  mutate(across(c(BAL_EOS_no_V5, BAL_LYM_no_V5, BAL_PMN_no_V5, 
                  BAL_MONO_no_V5), ~./BAL_volume.ml_V5/1E3)) %>% 
  select(-BAL_volume.ml_V4, -BAL_volume.ml_V5) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("group", "variable", "unit", "visit"), sep="_") %>% 
  select(-unit) %>% 
  filter(variable %in% c("WBC","EOS")) %>% 
  mutate(variable = recode_factor(factor(variable),
                                  "WBC"="WBC total",
                                  "EOS"="EOS per ul")) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_point() +
  geom_line(aes(group=ptID), color="grey") +
  facet_wrap(group~variable, scales="free", ncol=2) +
  theme_classic(base_size = 20) +
  labs(x="", y="")

ggsave("figs/for_ppt/cells.png", plot1, width=8, height=8)

#### Asthma ####
asthma <- read_csv("data_raw/addtl.data/P337_patient.metadata.csv")

plot2 <- asthma %>% 
  select(ptID, FeNO.PreBro_V4, FeNO.PreBro_V5) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("variable","visit"), sep="_") %>% 
  select(ptID, visit, variable, value) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_point() +
  geom_line(aes(group=ptID), color="grey") +
  facet_wrap(~variable, scales="free", ncol=2) +
  theme_classic(base_size = 20) +
  labs(x="", y="")

ggsave("figs/for_ppt/FeNO.png", plot2, width=4, height=8)

plot3 <- asthma %>% 
  select(ptID, FEV1.pctPP.preAlbuterol_V4:FEV1.pct.rev_V5) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("variable","visit"), sep="_") %>% 
  select(ptID, visit, variable, value) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_point() +
  geom_line(aes(group=ptID), color="grey") +
  facet_wrap(~variable, scales="free", ncol=2) +
  theme_classic(base_size = 20) +
  labs(x="", y="")

ggsave("figs/for_ppt/FEV1.png", plot3, width=8, height=4)

#### Cytokines ####
elisa <- read_csv("data_raw/addtl.data/P337_BAL.ELISA.csv")

plot4 <-elisa %>% 
  filter(stim =="stim") %>% 
  select(-stim,-unit) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("variable","visit"), sep="_") %>% 
  select(ptID, visit, variable, value) %>% 
  #Convert to log10
  mutate(value = ifelse(value==0,0, log10(value))) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_point() +
  geom_line(aes(group=ptID), color="grey") +
  facet_wrap(~variable, scales="free", nrow=2) +
  theme_classic(base_size = 20) +
  labs(y="log10( pg/ml )", x="")

ggsave("figs/for_ppt/elisa.png", plot4, width=16, height=8)

#### multiplex ####
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv")

plot5 <- plex %>% 
  #Correct error
  rename(FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("variable","visit"), sep="_") %>% 
  select(ptID, visit, variable, value)  %>% 
  #Convert to log10
  mutate(value = ifelse(value==0,0, log10(value))) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_point() +
  geom_line(aes(group=ptID), color="grey") +
  facet_wrap(~variable, scales="free_y", ncol=14) +
  theme_classic(base_size = 12) +
  labs(y="log10( pg/ml )", x="")

#plot5
ggsave("figs/for_ppt/mplex.png", plot5, width=20, height=8)

#### fMRI ####
library(readxl)

#Time point 1 = visit 4
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

plot6 <- neuro %>% 
  filter(neuro %notin% c("LSI","BDI")) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_point() +
  geom_line(aes(group=donorID), color="grey") +
  facet_wrap(~neuro, scales="free", nrow=2) +
  theme_classic(base_size = 20) +
  labs(x="", y="")
#plot6
ggsave("figs/for_ppt/fMRI.png", plot6, width=12, height=6)


#### Correlation ####
library(corrplot)
R <- read_csv("results/correlation/R_pearson_P337_delta.csv")
P <- read_csv("results/correlation/P_pearson_P337_delta.csv")

#List signif variables
modules <- colnames(R)[grepl("BAL_EOS.|BAL_PMN.", colnames(R))]
pct <- c("EOS.pct", "PMN.pct")
fMRI <- colnames(R)[grepl("insula|HO|ACC", colnames(R))]
multi <- colnames(R)[grepl("multiplex", colnames(R))]

all_vars <- c(modules,pct,fMRI,multi)

#### Calculate FDR ####
P.lower <-  P %>% 
  column_to_rownames("variable")
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

#### Signif corr ####
#Signig for EOS mod 1 or 2
multi.signif <- FDR %>% 
  filter(variable %in% c("BAL_EOS.pct_01","BAL_EOS.pct_02")) %>% 
  select(variable, all_of(multi)) %>% 
  pivot_longer(-variable) %>% 
  filter(value <= 0.05) %>% 
  arrange(variable, value) %>% 
  distinct(name) %>% unlist(use.names = FALSE)

#all_vars_signif <- c(modules,pct,fMRI,multi.signif[-20], "VEGF.multiplex")
all_vars_signif <- all_vars

R.signif <- R %>% 
  select(variable, all_of(all_vars_signif))
R.signif <- R.signif %>% 
  filter(variable %in% all_vars_signif) %>% 
  mutate(variable = factor(variable, levels=colnames(R.signif)[-1])) %>% 
  arrange(variable)

FDR.signif <- FDR %>% 
  select(variable, all_of(all_vars_signif))
FDR.signif <- FDR.signif %>% 
  filter(variable %in% all_vars_signif) %>% 
  mutate(variable = factor(variable, levels=colnames(FDR.signif)[-1])) %>% 
  arrange(variable)

# Check
identical(colnames(R.signif)[-1], as.character(R.signif$variable))
identical(R.signif$variable, FDR.signif$variable)
identical(colnames(R.signif), colnames(FDR.signif))

png(file = "figs/for_ppt/corr.fdr.png", 
    height=12, width=12, units="in", res=150)

corrplot(as.matrix(column_to_rownames(R.signif, "variable")),
         method="color", order="original", type="upper", diag=FALSE,
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90, cl.lim=c(-1,1),
         #Change colorlegend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = as.matrix(column_to_rownames(FDR.signif, "variable")),
         insig = "label_sig",
         sig.level = c(0.01,0.05), 
         pch.cex=0.9, pch.col="white",
         #Add title
         title="FDR *<0.05 **<0.01", 
         mar = c(0,0,2,0))

dev.off()

#### fMRI vs neuro ####
R.signif <- R %>% 
  filter(variable %in% c("BDI","LSI")) %>% 
  select(variable, R_vA_insula:FeNO.PreBro,-BDI, -LSI )

P.signif <- P %>% 
  filter(variable %in% c("BDI","LSI")) %>% 
  select(variable, R_vA_insula:FeNO.PreBro,-BDI, -LSI )

png(file = "figs/for_ppt/corr2.png", height=4, width=5, units="in", res=150)

corrplot(as.matrix(column_to_rownames(R.signif, "variable")),
         method="color", order="original",
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90, cl.lim=c(-1,1),
         #Change colorlegend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = as.matrix(column_to_rownames(P.signif, "variable")),
         insig = "label_sig",
         sig.level = c(0.01, 0.05), 
         pch.cex=0.9, pch.col="white",
         #Add title
         title="P *<0.05 **<0.01", 
         mar = c(0,0,2,0))

dev.off()

#### dotplot corr ####
to.keep.pair <- to.keep.var %>% 
  pivot_longer(-variable) %>% 
  filter(value <= 0.05)

delta <- read_csv("data_clean/P337.delta.all.csv")

#colors
col.vec <- c('#33a02c','#a6cee3','#b2df8a','#fdbf6f','#fb9a99',
             '#ff7f00','#e31a1c','#cab2d6','#6a3d9a','#1f78b4')
names(col.vec) <- unique(to.keep.pair$variable)

plot.ls.pos <- list()
plot.ls.neg <- list()

for(i in 1:nrow(to.keep.pair)){

  R.val <- R %>%
    filter(variable == to.keep.pair$variable[i]) %>%
    select(to.keep.pair$name[i]) %>%
    unlist(use.names = FALSE) %>%
    round(., digits=3)
  P.val <- P %>%
    filter(variable == to.keep.pair$variable[i]) %>%
    select(to.keep.pair$name[i]) %>%
    unlist(use.names = FALSE)%>%
    round(., digits=3)

  plot.temp <- delta %>%

    ggplot(aes_string(y=to.keep.pair$variable[i],
               x=to.keep.pair$name[i])) +
    geom_point(color=col.vec[to.keep.pair$variable[i]], size=2) +
    geom_smooth(method="lm", se=FALSE, color="darkgrey")+
    theme_classic() +
    labs(y=paste("Delta log2 expression\n", to.keep.pair$variable[i]),
         x=paste("Delta",to.keep.pair$name[i]),
         title=paste("R^2 =", R.val, ", P =", P.val)) +
    theme(plot.title = element_text(size=10))

  if(R.val < 0){
    plot.ls.neg[[as.character(i)]] <- plot.temp
  } else{
    plot.ls.pos[[as.character(i)]] <- plot.temp
  }
}

plot7a <- plot_grid(plotlist=plot.ls.pos)
plot7b <- plot_grid(plotlist=plot.ls.neg)

ggsave("figs/for_ppt/corr.dotplota.png", plot7a, width=12, height=8)
ggsave("figs/for_ppt/corr.dotplotb.png", plot7b, width=12, height=8)

#####
vars <- colnames(R)[c(3,9,11,15,17,18,21,24,25,32,38:47,95)]
R.val <- R %>%
  filter(variable %in% vars) %>%
  select(variable, all_of(vars)) 
P.val <- P %>%
  filter(variable %in% vars) %>%
  select(variable, all_of(vars)) 

  
corrplot(as.matrix(column_to_rownames(R.val, "variable")),
         method="color", order="original",
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90, cl.lim=c(-1,1),
         #Change colorlegend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = as.matrix(column_to_rownames(P.val, "variable")),
         insig = "label_sig",
         sig.level = c(0.01, 0.05), 
         pch.cex=0.9, pch.col="white",
         #Add title
         title="P *<0.05 **<0.01", 
         mar = c(0,0,2,0))

