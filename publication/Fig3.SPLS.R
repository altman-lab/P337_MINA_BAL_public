library(tidyverse)
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
                       Z=c(5.5,7.4,27.4,NA)) %>% 
  mutate(facet.name = fct_relevel(factor(facet.name), rev))

#Order labels and variables
##Both components
X.both <- c("CCL2","BAL EOS 02")
Y.both <- c("L vA insula","R HO amgy")
##Comp1 only
X1 <- loadX %>% 
  filter(comp1 != 0 & var %notin% X.both) %>% 
  arrange(-comp1) %>% 
  select(var) %>% unlist(use.names=FALSE)
Y1 <- loadY %>% 
  filter(comp1 != 0 & var %notin% Y.both) %>% 
  arrange(-comp1) %>% 
  select(var) %>% unlist(use.names=FALSE)
##Comp2 only
X2 <- loadX %>% 
  filter(comp2 != 0 & var %notin% X.both) %>% 
  arrange(-comp2) %>% 
  select(var) %>% unlist(use.names=FALSE)
Y2 <- loadY %>% 
  filter(comp2 != 0 & var %notin% Y.both) %>% 
  arrange(-comp2) %>% 
  select(var) %>% unlist(use.names=FALSE)

plot.dat <- loadXY %>% 
  mutate(comp = recode(comp, "comp1"="Component 1", "comp2"="Component 2")) %>% 
  mutate(facet.name = fct_relevel(factor(facet.name), rev)) %>% 
  #order variables
  mutate(var = factor(var, levels=c(X2,X1,X.both,Y2,Y1,Y.both))) %>%
  arrange(desc(var))

#### Plot ####
plot <- plot.dat %>% 
  ggplot(aes(x=var, y=value, fill=comp)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"),
           stat="identity") +
  scale_x_reordered() +
  coord_flip() + 
  theme_classic() +
  facet_wrap(~facet.name, scales="free_y", ncol=2) +
  labs(x="Selected variable", y="Loading", fill="") +
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("#7CAE00","#C77CFF")) +
  theme(legend.position = "bottom") +
  geom_vline(data=line.dat, aes(xintercept=Z), color="grey")
plot

#### Save ####
ggsave("publication/Fig3.SPLS.png", plot, width=6, height=5)

#### Heatmap ####
library(mixOmics)

cols <- colorRampPalette(c("#2166ac","white","#b2182b"))
cim(result.spls, margins = c(7, 7), 
    color=cols(20),
    save="png", name.save = "publication/Fig3.SPLS.heatmap")
