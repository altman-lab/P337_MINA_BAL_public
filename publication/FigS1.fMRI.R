library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)

#### Data ####
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
  full_join(neuro) %>% 
  filter(neuro != "BDI" & neuro != "LSI") %>% 
  mutate(neuro = gsub("_"," ",neuro)) %>%  
  mutate(visit = recode_factor(factor(visit), "V4"="Pre","V5"="Post"))

#### Heatmap delta ####
plot.dat <- neuro %>% 
  mutate(donorID = paste("MA",idnum, sep="")) %>% 
  pivot_wider(names_from = visit) %>% 
  mutate(delta=Post-Pre) %>% 
  #remove fxnal metrics
  filter(!grepl("3way",neuro) & !grepl("LPR", neuro)) %>% 
  #fix name
  mutate(neuro = gsub("amgy","amyg",neuro)) %>% 
  #wide format
  dplyr::select(-idnum,-Pre,-Post) %>% 
  pivot_wider(names_from = neuro, values_from = delta) %>% 
  column_to_rownames("donorID") %>% 
  t()

#Legend
col_fun<- colorRamp2(c(-100, 0, 100), c("blue", "white", "red"))
lgd <- Legend(col_fun = col_fun, title = "Post - Pre challenge",
              direction = "horizontal", title_position = "topcenter")

plot <- Heatmap(plot.dat, show_heatmap_legend = FALSE)

png("publication/FigS1.fMRI.png", height=3.5, width=6, units="in", res=150)
  draw(plot, padding=unit(c(0,0,15,0),"mm"))
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.99, "npc"), just = c("center", "top"))
dev.off()


pdf("publication/FigS1.fMRI.pdf", height=3.5, width=6)
draw(plot, padding=unit(c(0,0,15,0),"mm"))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.99, "npc"), just = c("center", "top"))
dev.off()
