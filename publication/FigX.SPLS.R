library(tidyverse)
library(drlib)
library(ggrepel)
library(cowplot)
set.seed(3698)

#### Data ####
load("results/PLS/SPLS.RData")

#### Loadings ####
loadX <- as.data.frame(result.spls$loadings$X) %>% 
  rownames_to_column("var") %>% 
  mutate(space="X")

loadXY <- as.data.frame(result.spls$loadings$Y) %>% 
  rownames_to_column("var") %>% 
  mutate(space="Y") %>% 
  bind_rows(loadX) %>% 
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
  
#### Loading plot1 ####
plot1 <- loadXY %>% 
  #reorder facets
  mutate(name.index = as.character(as.numeric(as.factor(name))),
         name.index = factor(name.index, levels=c(3,4,1,2)),
         name.index = as.numeric(as.factor(name.index))) %>% 
  mutate(name = fct_reorder(name, name.index)) %>% 
  
  ggplot(aes(x=reorder_within(var, -value, list(name, space)), 
             y=value, fill=col.group)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"),
           stat="identity") +
  scale_x_reordered() +
  coord_flip() + 
  theme_classic() +
  facet_wrap(~name, scales="free_y", ncol=2) +
  labs(x="Selected variable", y="Loading", fill="") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0)
plot1

#### Loading plot2 ####
plot2 <- loadXY %>% 
  ggplot(aes(x=reorder_within(var, -value, list(space)), 
             y=value, fill=name)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"),
           stat="identity") +
  scale_x_reordered() +
  coord_flip() + 
  theme_classic() +
  facet_wrap(~space, scales="free_y", ncol=2) +
  labs(x="Selected variable", y="Loading", fill="") +
  theme(legend.key.size = unit(0.8, "cm"))
plot2

#### Loading plot3 ####

plot3 <- loadXY %>% 
  select(-name) %>% 
  pivot_wider(names_from = comp) %>% 
  replace_na(list(comp1=0,comp2=0)) %>% 
  
  ggplot(aes(x=comp1, y=comp2)) +
  geom_point(aes(color=col.group), size=2) +
  geom_text_repel(aes(label=var), direction="both",
                  min.segment.length = unit(0, 'lines'),
                  show.legend = FALSE, size=3, max.overlaps = 100) +
  coord_fixed() + 
  theme_classic() +
  labs(x=paste("Component 1\nmodule + cytokine = ", 
               round(result.spls$explained_variance$X[1]*100, digits=1), 
               "%\nfMRI = ",
               round(result.spls$explained_variance$Y[1]*100, digits=1),
               "%", sep=""), 
       y=paste("Component 2\nmodule + cytokine = ", 
               round(result.spls$explained_variance$X[2]*100, digits=1),
               "%\nfMRI = ",
               round(result.spls$explained_variance$Y[2]*100, digits=1),
               "%", sep=""),
       color="Dataset") +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed")

plot3

#### Save ####

ggsave("publication/FigX.SPLS1.png", plot1, width=8, height=7)
ggsave("publication/FigX.SPLS2.png", plot2, width=8, height=6)
ggsave("publication/FigX.SPLS3.png", plot3, width=6, height=7)
