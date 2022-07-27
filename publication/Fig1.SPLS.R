library(tidyverse)
library(ggforce)
library(readxl)
library(Hmisc)
library(cowplot)
set.seed(3698)
`%notin%` <- Negate(`%in%`)

#### PLS Data ####
load("results/PLS/SPLS.RData")

#### Loadings ####
loadX <- as.data.frame(result.spls$loadings$X) %>% 
  rownames_to_column("var") %>% 
  mutate(space="X")
loadY <- as.data.frame(result.spls$loadings$Y) %>% 
  rownames_to_column("var")%>% 
  mutate(space="Y")

loadXY <- bind_rows(loadY, loadX) %>% 
  pivot_longer(comp1:comp2, names_to = "comp") %>% 
  mutate(name = ifelse(space=="X" & comp=="comp1", 
                       paste("modules + cytokines\nComponent 1 = ", 
                             round(result.spls$explained_variance$X[1]*100, digits=1), 
                             "%", sep=""),
                       ifelse(space=="X" & comp=="comp2", 
                              paste("modules + cytokines\nComponent 2 = ",
                                    round(result.spls$explained_variance$X[2]*100, digits=1),
                                    "%", sep=""),
                              ifelse(space=="Y" & comp=="comp1", 
                                     paste("fMRI\nComponent 1 = ",
                                           round(result.spls$explained_variance$Y[1]*100, digits=1), 
                                           "%", sep=""),
                                     ifelse(space=="Y" & comp=="comp2", 
                                            paste("fMRI\nComponent 2 = ",
                                                  round(result.spls$explained_variance$Y[2]*100, digits=1),
                                                  "%", sep=""),
                                            NA))))) %>% 
  mutate(facet.name = ifelse(space=="X", "modules + cytokines",
                             ifelse(space=="Y", 
                                    paste("fMRI"),
                                    NA))) %>% 
  filter(value != 0) %>% 
  #color group
  mutate(col.group = ifelse(grepl("BAL EOS", var), "module",
                            ifelse(grepl("L |R |pACC|dACC", var), "fMRI",
                                   "cytokine"))) %>% 
  mutate(var = ifelse(var=="BAL EOS 02", "EOS02 module", var)) %>% 
  #fix protein names
  mutate(var = gsub("^IL","IL-",var),
         var = recode(var,
                      "FLT3L"="FLT3LG",
                      "TGFA"="TGF-A",
                      "GCSF"="G-CSF",
                      "PDGFAA"="PDGF-AA",
                      "IFNA2"="IFN-A2",
                      "GMCSF"="GM-CSF"))


#### Plot setup ####
#cut lines for components mapped to
line.dat <- data.frame(facet.name = c("modules + cytokines","fMRI"), 
                       Z=c(2.5,5.5)) %>% 
  mutate(facet.name = fct_relevel(factor(facet.name), rev))

#Order labels and variables
##Comp1 only
X1 <- loadXY %>% 
  filter(space == "X" & comp == "comp1") %>% 
  arrange(-value) %>% 
  mutate(var=ifelse(var=="BAL EOS 02", "EOS02 module", var)) %>% 
  dplyr::select(var) %>% unlist(use.names=FALSE)
##Comp2 only
X2 <- loadXY %>% 
  filter(space == "X" & comp == "comp2") %>% 
  arrange(-value) %>%
  dplyr::select(var) %>% unlist(use.names=FALSE)

plot.dat <- loadXY %>% 
  mutate(comp = recode_factor(factor(comp), "comp2"="Component 2", 
                              "comp1"="Component 1")) %>% 
  mutate(facet.name = fct_relevel(factor(facet.name), rev)) %>% 
  #order variables
  mutate(var = factor(var, levels=unique(c(X2,X1,
                                           "R vA insula","R dA insula",
                                           "dACC","pACC","R HO amyg",
                                           "L dA insula","L HO amyg",
                                           "L vA insula")))) %>%
  arrange(desc(var)) 

#### Plot PLS ####
plot <- plot.dat %>% 
  ggplot(aes(x=var, y=value, fill=comp)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"),
           stat="identity") +
  #scale_x_reordered() +
  coord_flip() + 
  theme_classic() +
  ggforce::facet_col(facets = vars(facet.name), 
                     scales = "free_y", 
                     space = "free") +
  labs(x="Selected variable", y="Loading", fill="") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("#7CAE00","#C77CFF")) +
  theme(legend.position = "bottom") +
  geom_vline(data=line.dat, aes(xintercept=Z), color="grey") +
  guides(fill = guide_legend(reverse=TRUE))
plot

#### Heatmap Data ####
#### fMRI data ####
#Time point 1 = visit 4
neuro <- read_excel(sheet="T1",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(visit="V4")

#Time point 2 = visit 5
neuro <- read_excel(sheet="T2",
                    "data_raw/addtl.data/extraced.clusters.Matt.Altman_wbaseline_psychdata.xlsx") %>% 
  #long format
  pivot_longer(-idnum, names_to="neuro") %>% 
  #add sample variable
  mutate(visit="V5") %>% 
  #Combine with other visit
  full_join(neuro)  %>% 
  #Format idnum to match RNAseq data
  mutate(idnum = paste("MA",idnum, sep="")) %>% 
  rename(donorID=idnum) %>% 
  #remove non fMRI vars
  filter(neuro != "BDI" & neuro != "LSI") %>% 
  #remove fxnal metrics
  filter(!grepl("3way",neuro) & !grepl("LPR", neuro)) %>% 
  #Remove failed sample
  filter(donorID != "MA1012") %>% 
  #Clean names
  mutate(neuro = gsub("_"," ",neuro)) %>%  
  mutate(visit = recode_factor(factor(visit), "V4"="Pre","V5"="Post")) %>% 
  mutate(neuro = gsub("amgy","amyg",neuro))

# Calculate delta
neuro.delta <- neuro %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=Post-Pre)

#### Cell data ####
load("data_clean/P337_BAL_data.RData")

cell.delta <- dat.BAL.abund.norm.voom$targets %>% 
  select(visit, donorID, EOS.pct, PMN.pct) %>% 
  pivot_longer(EOS.pct:PMN.pct) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=V5-V4) %>% 
  select(-V4, -V5) %>% 
  pivot_wider(values_from = delta)

##### Module data ####
attach("data_clean/P337_BAL_module_data.RData")

mod.delta <- mod.voom %>% 
  pivot_longer(-module, names_to = "libID") %>% 
  left_join(select(dat.BAL.abund.norm.voom$targets, libID, visit, donorID)) %>% 
  select(-libID) %>% 
  filter(grepl("EOS.pct_02", module)) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(EOS02=V5-V4) %>% 
  select(-V4, -V5, -module)

##### protein data ####
#SPLS-selected
attach("results/PLS/SPLS.RData")

cyto <- as.data.frame(result.spls$loadings$X) %>% 
  filter(comp1 !=0 | comp2 != 0)

plex.delta<- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  pivot_longer(-ptID) %>% 
  separate(name, into=c("name","visit")) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=V5-V4) %>% 
  select(-V4, -V5) %>% 
  filter(name %in% rownames(cyto)) %>% 
  #fix protein names
  mutate(name = gsub("^IL","IL-",name),
         name = recode(name,
                      "FLT3L"="FLT3LG",
                      "TGFA"="TGF-A",
                      "GCSF"="G-CSF",
                      "PDGFAA"="PDGF-AA",
                      "IFNA2"="IFN-A2",
                      "GMCSF"="GM-CSF"))%>% 
  pivot_wider(values_from = delta) %>% 
  rename(donorID=ptID)

#### combine data ####
X <- neuro.delta %>% 
  select(donorID, neuro, delta) %>% 
  pivot_wider(names_from = neuro, values_from = delta) %>% 
  full_join(cell.delta) %>% 
  full_join(mod.delta) %>% 
  full_join(plex.delta) %>% 
  rename(`BAL EOS 02`=EOS02, `EOS %`=EOS.pct, `PMN %`=PMN.pct) %>% 
  rename_all(~gsub("_"," ",.)) %>% 
  #order Y
  select(donorID, #`EOS %`, `PMN %`,
         `FLT3LG`,`IL-17A`,`BAL EOS 02`,`G-CSF`,`IL-10`,
         `PDGF-AA`,`CCL27`,`IFN-A2`,`GM-CSF`,`CCL21`,
         `TGF-A`,`CCL26`,`IL-23`,`CCL7`,`CXCL1`,`IL-16`,
         `L vA insula`,`L HO amyg`,`L dA insula`,`R HO amyg`,
         `pACC`,`dACC`,`R dA insula`,`R vA insula`) %>% 
  column_to_rownames("donorID") %>% 
  as.matrix()

#### Correlation ####
cor.result <- Hmisc::rcorr(X, type = "pearson")

#Arrange as in barplot
cor.arrange <- as.data.frame(cor.result$r) %>% 
  rownames_to_column("X_variable") %>% 
  pivot_longer(-X_variable, values_to = "r")
cor.arrange <- as.data.frame(cor.result$P) %>% 
  rownames_to_column("X_variable") %>% 
  pivot_longer(-X_variable, values_to = "p") %>% 
  full_join(cor.arrange) %>% 
  #order Y
  mutate(name = factor(name, levels=rev(colnames(X))),
         X_variable = factor(X_variable, 
                             levels=c(colnames(X)))) %>% 
  rowwise() %>%
  mutate(pair = sort(c(X_variable, name)) %>% paste(collapse = ",")) %>%
  group_by(pair) %>%
  distinct(pair, .keep_all = T) %>% 
  filter(X_variable != name) %>% 
  ungroup() %>% 
  #signif label
  mutate(label = ifelse(p < 0.05, "*",""),
         FDR = p.adjust(p),
         label2 = ifelse(FDR < 0.05, "*",""))

#### plot ####
plot2 <- cor.arrange %>% 
  ggplot(aes(x=X_variable, y=name, fill=r)) +
  geom_tile() +
  geom_text(aes(label = label)) +
  labs(x="", fill="Pearson") +
  # coord_flip() +
  coord_fixed()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.position = c(0.7,0.8), legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5,
                               title.position = "top", title.hjust = 0.5)) +
  scale_fill_gradient2(low="darkblue",mid="white",high="darkred",
                       limits=c(-1,1))
plot2

#### Save ####
ggsave("publication/Fig1.SPLS.png", 
       plot_grid(plot,plot2, rel_widths = c(0.5, 1), 
                 labels = c("(A)","(B)"),
                 align="hv", axis="tb", nrow = 1),
       width=8.5, height=5)
ggsave("publication/Fig1.SPLS.pdf", 
       plot_grid(plot,plot2, rel_widths = c(0.5, 1), 
                 labels = c("(A)","(B)"),
                 align="hv", axis="tb", nrow = 1),
       width=8.5, height=5)
