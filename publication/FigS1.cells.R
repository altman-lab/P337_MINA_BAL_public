library(tidyverse)

#### Data ####
load("data_clean/P337_BAL_data.RData")

#### Stats ####
cell.test <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, ends_with(".pct")) %>% 
  pivot_longer(-c(donorID, visit)) %>% 
  arrange(name, donorID, visit) %>% 
  group_by(name, visit) %>% 
  summarise(value = list(value)) %>% 
  pivot_wider(names_from = visit) %>% 
  group_by(name) %>% 
  mutate(p = t.test(unlist(V4), unlist(V5))$p.value,
         t = t.test(unlist(V4), unlist(V5))$statistic) %>% 
  ungroup() %>% 
  mutate(FDR = p.adjust(p))

#### plot ####
p1 <- dat.BAL.abund.norm.voom$targets %>% 
  select(donorID, visit, ends_with(".pct")) %>% 
  pivot_longer(-c(donorID, visit)) %>% 
  group_by(name, donorID) %>% 
  mutate(diff = value[visit=="V5"]-value[visit=="V4"],
         diff.col = ifelse(diff<0, "down","up"),
         diff.col = factor(diff.col, levels=c("up","down"))) %>%
  ungroup() %>% 
  left_join(cell.test) %>% 
  mutate(cell = gsub(".pct","",name),
         cell = recode(cell, 
                       "MONO"="monocyte",
                       "EOS"="eosinophil", 
                       "LYM"="lymphocyte",
                       "Epi"="epithelial"),
         FDR.format = ifelse(FDR > 0.01, round(FDR, digits=3),
                             formatC(FDR, digits = 2, format="e")),
         cell = paste(cell, paste("FDR =", FDR.format), sep="\n"),
         visit = recode_factor(factor(visit), "V4"="Pre", "V5"="Post")) %>% 
  
  ggplot(aes(x=visit, y=value)) +
  geom_path(aes(group=donorID, color=diff.col)) +
  geom_point() +
  theme_classic() +
  labs(x="Segmented bronchial provocation with allergen (SBP-Ag)",
       y="Cell %", color="Post - Pre\nchange") +
  facet_grid(~cell, scales="free_x") +
  scale_color_manual(values=c("down"="#74add1","up"="#d73027"))
p1

#### Save ####
ggsave("publication/FigS1.cells.png", p1, height=2.5, width=7.3)

ggsave("publication/FigS1.cells.pdf", p1, height=2.5, width=7.2)
