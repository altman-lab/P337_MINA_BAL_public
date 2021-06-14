library(tidyverse)
library(patchwork)
library(mixOmics)

#### Data ####
attach("results/PLS/SPLS.RData")
attach("results/PLS/spls.tune.RData")
attach("results/PLS/spls.MAEzoom.tune.RData")

#### axis loading ####
plot.load <- as.data.frame(spls$explained_variance) %>% 
  mutate(Component = c(1:10)) %>% 
  pivot_longer(X:Y, names_to = "space") %>% 
  #mutate(space = recode(space, "X"="X: gene modules + proteins","Y"="Y: fMRI")) %>% 
  
  ggplot(aes(x=Component, y=value)) +
  geom_point(size=2) +
  facet_wrap(~space) +
  theme_classic() +
  labs(y="Loading") +
  scale_x_continuous(breaks=seq(1,10,1)) +
  geom_vline(xintercept = 2, lwd=1, lty="dashed")

#### MAE ####
plot <- plot(spls.zoom.MAE.X, legend.position = 'topright') +
  labs(y="Mean absolute error\n(MAE)", color="Component",
       title="X: gene modules + proteins") +
  geom_vline(xintercept = 2, lwd=1, lty="dashed", color="#FF6633") +
  geom_vline(xintercept = 14, lwd=1, lty="dashed", color="#3366CC") +
  theme(legend.position = "bottom")
#plot

plot2 <- plot(spls.MAE.Y, legend.position = 'topright') +
  labs(y="Mean absolute error\n(MAE)", color="Component",
       title="Y: fMRI") +
  geom_vline(xintercept = 8, lwd=1, lty="solid", color="#3366CC") +
  geom_vline(xintercept = 8, lwd=1, lty="dashed",color="#FF6633") +
  theme(legend.position = "none")
#plot2

#### Save ####
layout <- "
ABB
CCC
"

plot.all <- plot.load+plot2+plot +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = '(',tag_suffix = ')') + 
  plot_layout(design=layout)
ggsave(plot.all, filename = "publication/FigSX.MAE.png", height=6, width=8)
ggsave(plot.all, filename = "publication/FigSX.MAE.pdf", height=6, width=8)
