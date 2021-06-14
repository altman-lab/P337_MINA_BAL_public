library(tidyverse)
library(readxl)
library(drlib)
library(ggrepel)
library(cowplot)
set.seed(3698)
`%notin%` <- Negate(`%in%`)

#### Data ####
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
  mutate(facet.name = ifelse(space=="X", 
                       paste("modules + cytokines\nComponent 1 = ", 
                             round(result.spls$explained_variance$X[1]*100, digits=1), 
                             "%\nComponent 2 = ",
                             round(result.spls$explained_variance$X[2]*100, digits=1),
                             "%", sep=""),
                       ifelse(space=="Y", 
                              paste("fMRI\nComponent 1 = ",
                                    round(result.spls$explained_variance$Y[1]*100, digits=1), 
                                    "%\nComponent 2 = ",
                                    round(result.spls$explained_variance$Y[2]*100, digits=1),
                                    "%", sep=""),
                              NA))) %>% 
  filter(value != 0) %>% 
  #color group
  mutate(col.group = ifelse(grepl("BAL EOS", var), "module",
                            ifelse(grepl("L |R |pACC|dACC", var), "fMRI",
                                   "cytokine"))) 
  
#### Plot setup ####
#cut lines for components mapped to
line.dat <- data.frame(facet.name = c(paste("modules + cytokines\nComponent 1 = ", 
                                   round(result.spls$explained_variance$X[1]*100,
                                         digits=1), 
                                   "%\nComponent 2 = ",
                                   round(result.spls$explained_variance$X[2]*100,
                                         digits=1),
                                   "%", sep=""),
                             paste("fMRI\nComponent 1 = ",
                                   round(result.spls$explained_variance$Y[1]*100,
                                         digits=1), 
                                   "%\nComponent 2 = ",
                                   round(result.spls$explained_variance$Y[2]*100,
                                         digits=1),
                                   "%", sep="")), 
                       Z=c(2.5,5.5)) %>% 
  mutate(facet.name = fct_relevel(factor(facet.name), rev))

#Order labels and variables
##Comp1 only
X1 <- loadX %>% 
  filter(comp1 != 0 ) %>% 
  arrange(-comp1) %>% 
  dplyr::select(var) %>% unlist(use.names=FALSE)
##Comp2 only
X2 <- loadX %>% 
  filter(comp2 != 0) %>% 
  arrange(-comp2) %>% 
  dplyr::select(var) %>% unlist(use.names=FALSE)

plot.dat <- loadXY %>% 
  mutate(comp = recode_factor(factor(comp), "comp2"="Component 2", "comp1"="Component 1")) %>% 
  mutate(facet.name = fct_relevel(factor(facet.name), rev)) %>% 
  #order variables
  mutate(var = factor(var, levels=unique(c(X2,X1,
                                           "R vA insula","R dA insula","dACC","pACC","R HO amyg",
                                           "L dA insula","L HO amyg","L vA insula")))) %>%
  arrange(desc(var))

#### Plot ####
plot <- plot.dat %>% 
  ggplot(aes(x=var, y=value, fill=comp)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"),
           stat="identity") +
  #scale_x_reordered() +
  coord_flip() + 
  theme_classic() +
  facet_wrap(~facet.name, scales="free_y", ncol=2) +
  labs(x="Selected variable", y="Loading", fill="") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("#7CAE00","#C77CFF")) +
  theme(legend.position = "bottom") +
  geom_vline(data=line.dat, aes(xintercept=Z), color="grey") +
  guides(fill = guide_legend(reverse=TRUE))
plot

#### Heatmap ####
#Calculate post-pre change (delta)
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
  dplyr::select(donorID, visit, everything(), -idnum) %>% 
  arrange(donorID, visit) 

#fMRI
neuro.delta <- neuro %>% 
  #transpose
  pivot_longer(-c(donorID,visit)) %>% 
  pivot_wider(names_from = visit) %>% 
  #calculate delta
  rowwise() %>% 
  mutate(delta = post-pre) %>% 
  dplyr::select(-pre,-post) %>% 
  pivot_wider(values_from = delta) 

#modules
attach("data_clean/P337_BAL_data.RData")
attach("data_clean/P337_BAL_module_data.RData")

mod.delta <- mod.voom %>% 
  #Remove modules 00
  filter(!grepl("_00",module)) %>% 
  #Add metadata to denote V4, V5
  pivot_longer(-module, names_to = "libID") %>% 
  full_join(dplyr::select(dat.BAL.abund.norm.voom$targets, 
                          libID, donorID, visit), by = "libID") %>% 
  #separate V4 and V5
  dplyr::select(-libID) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  #shorten module names
  mutate(module=gsub("P337_", "", module),
         module=gsub(".pct", "", module)) %>% 
  #wide format
  dplyr::select(-V4,-V5) %>% 
  arrange(donorID, module) %>% 
  pivot_wider(names_from = module, values_from = delta)

#Protein
plex.delta <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("name","visit"), sep="_") %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta = V5-V4) %>% 
  drop_na(delta) %>% 
  dplyr::select(-V4,-V5) %>% 
  arrange(name)  %>% 
  pivot_wider(values_from = delta) 

#Combine X vars
X <- full_join(mod.delta, plex.delta) %>% 
  filter(donorID %in% neuro.delta$donorID) %>% 
  arrange(donorID) %>% 
  #matrix format
  column_to_rownames("donorID") %>% 
  as.matrix()

#Subset Y to donors in X
Y <- neuro.delta %>% 
  filter(donorID %in% rownames(X)) %>% 
  arrange(donorID) %>% 
  #matrix format
  column_to_rownames("donorID") %>% 
  as.matrix()

#Check
identical(rownames(X),rownames(Y))

cor <- cor(X,Y, method="pearson", use="pairwise.complete.obs") %>% 
  as.data.frame() %>% 
  rownames_to_column("X_variable")

#save
write_csv(cor, "results/PLS/pearson.correlation.csv")

#Arrange as in barplot
cor.arrange <- cor %>% 
  pivot_longer(-X_variable) %>% 
  mutate(X_variable = gsub("module_","",X_variable),
         X_variable = gsub("_"," ",X_variable),
         name = gsub("_"," ",name)) %>% 
  #Spls selected
  filter(X_variable %in% plot.dat$var & name %in% plot.dat$var) %>% 
  #order Y
  mutate(name = factor(name, levels=c(unique(as.character(plot.dat$var))[1:8]))) %>% 
  mutate(X_variable = factor(X_variable, 
                             levels=rev(unique(as.character(plot.dat$var))[9:24]))) %>% 
  arrange(X_variable, name) 

plot2 <- ggplot(cor.arrange, aes(x=X_variable, y=name, fill=value)) +
  geom_tile() +
  labs(x="", fill="Pearson") +
  scale_y_discrete(position = "right")  +
  coord_flip() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust=0),
    legend.position = "bottom",
    plot.margin = margin(0, 1, 0.2, 0.2, "cm"))+
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5,
                               title.position = "top", title.hjust = 0.5)) +
  scale_fill_gradient2(low="darkblue",mid="white",high="darkred",
                       limits=c(-0.8,0.8))
#plot2

#### Save ####
ggsave("publication/Fig1.SPLS.png", 
       plot_grid(plot,plot2, rel_widths = c(2.5,1), labels = c("(A)","(B)"),
                 align="hv", axis="b"),
       width=8, height=4)
ggsave("publication/Fig1.SPLS.pdf", 
       plot_grid(plot,plot2, rel_widths = c(2.5,1), labels = c("(A)","(B)"),
                 align="hv", axis="b"),
       width=8, height=4)
