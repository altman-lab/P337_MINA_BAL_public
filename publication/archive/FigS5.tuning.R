library(tidyverse)
library(patchwork)
library(mixOmics)

#### Data ####
attach("results/PLS/SPLS.RData")
attach("results/PLS/spls.tune2.RData")
attach("results/PLS/spls.tune2.zoom.RData")

colors <- c("#3366CC","#FF6633")

#loadings for selected # of features
load.select <-  data.frame(Component = c(1,2,1,2),
                           space = c("modules + cytokines","modules + cytokines","fMRI","fMRI"),
                           value = c(result.spls$prop_expl_var$X[1], 
                                     result.spls$prop_expl_var$X[2],
                                     result.spls$prop_expl_var$Y[1],
                                     result.spls$prop_expl_var$Y[2]))

#### axis loading ####
plot.load <- as.data.frame(spls$prop_expl_var) %>% 
  rownames_to_column() %>% 
  mutate(Component = as.numeric(gsub("comp","",rowname))) %>% 
  pivot_longer(X:Y, names_to = "space") %>% 
  mutate(space = recode(space, "X"="modules + cytokines", "Y"="fMRI")) %>% 

  ggplot(aes(x=Component, y=value)) +
  geom_point(size=1.5) +
  facet_wrap(~space) +
  theme_classic() +
  labs(y="Loading") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  #add actual selected loadings
  geom_point(data=load.select, aes(color=as.character(Component)), 
             size=2, shape=18, show.legend = FALSE) +
  scale_color_manual(values=colors)

# plot.load

#### MAE ####
plot2 <- spls.zoom.RSS.X$measure.pred %>% 
  filter(keepX <= 30) %>% 
  bind_rows(filter(spls.RSS.X$measure.pred, keepX>30)) %>% 
  filter(measure == "RSS" & V == "t" & comp < 3) %>% 
  ggplot() +
  aes(x = keepX, y=mean, color=as.character(comp)) +
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) + 
  geom_line(show.legend = FALSE) +
  labs(y="Residual Sum of Squares (RSS)", x="Number of selected features", 
       color="Component", title="modules + cytokines") +
  geom_vline(xintercept = 14, lwd=1, lty="dashed", color=colors[1]) +
  geom_vline(xintercept = 2, lwd=1, lty="dashed", color=colors[2]) +
  theme_classic() +
  scale_color_manual(values=colors)
plot2

plot4 <- spls.RSS.Y$measure.pred %>% 
  filter(measure == "RSS" & V == "t" & comp < 3) %>% 
  ggplot() +
  aes(x = keepX, y=mean, color=as.character(comp)) +
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) + 
  geom_line(show.legend = FALSE) +
  labs(y="Residual Sum of Squares (RSS)", x="Number of selected features", 
       color="Component", title="fMRI") +
  geom_vline(xintercept = 8, lwd=1,  color=colors[1]) +
  geom_vline(xintercept = 8, lwd=1, lty="dashed", color=colors[2]) +
  theme_classic() +
  scale_color_manual(values=colors)
# plot4

#### Save ####
layout <- "
AAABBB
CCCCCC
"

plot.all <- plot.load+plot4+plot2 +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = '(',tag_suffix = ')') + 
  plot_layout(design=layout)
ggsave(plot.all, filename = "publication/FigS5.tuning.png", height=6, width=8)
ggsave(plot.all, filename = "publication/FigS5.tuning.pdf", height=6, width=8)
