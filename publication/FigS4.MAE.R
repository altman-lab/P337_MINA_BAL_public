library(tidyverse)
library(patchwork)
library(mixOmics)

#### Data ####
attach("results/PLS/SPLS.RData")
attach("results/PLS/spls.tune.RData")
attach("results/PLS/spls.MAEzoom.tune.RData")

#loadings for selected # of features
load.select <-  data.frame(Component = c(1,2,1,2),
                           space = c("modules + cytokines","modules + cytokines","fMRI","fMRI"),
                           value = c(result.spls$explained_variance$X[1], 
                                     result.spls$explained_variance$X[2],
                                     result.spls$explained_variance$Y[1],
                                     result.spls$explained_variance$Y[2]))

#### axis loading ####
plot.load <- as.data.frame(spls$explained_variance) %>% 
  mutate(Component = c(1:10)) %>% 
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
  scale_color_manual(values=c("#3366CC","#FF6633"))


plot.load

#### MAE ####
plot <- plot(spls.zoom.MAE.X, legend.position = 'topright') +
  labs(y="Mean absolute error\n(MAE)", color="Component",
       title="modules + cytokines") +
  geom_vline(xintercept = 2, lwd=1, lty="dashed", color="#FF6633") +
  geom_vline(xintercept = 14, lwd=1, lty="dashed", color="#3366CC") +
  theme(legend.position = "bottom")
#plot

plot2 <- plot(spls.MAE.Y, legend.position = 'topright') +
  labs(y="Mean absolute error\n(MAE)", color="Component",
       title="fMRI") +
  geom_vline(xintercept = 8, lwd=1, lty="solid", color="#3366CC") +
  geom_vline(xintercept = 8, lwd=1, lty="dashed",color="#FF6633") +
  theme(legend.position = "none")
#plot2

#### Save ####
layout <- "
AAABBBB
CCCCCCC
"

plot.all <- plot.load+plot2+plot +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = '(',tag_suffix = ')') + 
  plot_layout(design=layout)
ggsave(plot.all, filename = "publication/FigS4.MAE.png", height=6, width=8)
ggsave(plot.all, filename = "publication/FigS4.MAE.pdf", height=6, width=8)
