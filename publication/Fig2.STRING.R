library(tidyverse)
`%notin%` <- Negate(`%in%`)
set.seed(3698)

#### STRING database ####
#STRING database
require(STRINGdb)
#Working with network objects
require(igraph) #vertex_attr()
#Graphing networks
require(ggraph)
require(scatterpie) #geom_scatterpie()
require(ggnetwork) #theme_blank()
require(scales)

string_db <- STRINGdb$new(version="11", species=9606,
                          score_threshold=700, input_directory="")

#### Genes in modules ####
load("data_clean/P337_BAL_module_data.RData")

mod.genes.temp <- mod.genes %>% 
    filter(grepl("EOS.pct_02", module)) %>% 
    distinct(hgnc_symbol) %>% unlist(use.names=FALSE)
  
#### Protein cytokines ####
load("results/PLS/SPLS.RData")

cyto <- as.data.frame(result.spls$loadings$X) %>% 
  rownames_to_column("var") %>% 
  filter(comp1 != 0 | comp2 != 0) %>% 
  filter(!grepl("BAL EOS", var)) %>% 
  mutate(gene = recode(gsub(".multiplex","",var),
                       "FLT3L"="FLT3LG",
                       "GCSF"="CSF3",
                       "IL12p70"="IL12A_IL12B",
                       "IL23"="IL23A",
                       "IP10"="CXCL10",
                       "PDGFAA"="PDGFA",
                       "PDGFAB"="PDGFA_PDGFB",
                       "sCD40L"="CD40LG")) %>% 
  #Unnest multi-annotations
  separate(gene, into=c("a","b"), sep="_", fill="right") %>% 
  pivot_longer(a:b, values_to = "gene") %>% 
  drop_na(gene)

#### Colors ####
# Set color values 
color.vec <- c("#fc8d62","#8da0cb")

dat.all <- data.frame(gene = c(mod.genes.temp, cyto$gene)) %>% 
  mutate(mod = ifelse(gene %in% mod.genes.temp, 1, 0),
         cyto = ifelse(gene %in% cyto$gene, 1, 0)) %>% 
  rowwise() %>% 
  mutate(total = sum(mod,cyto)) %>% 
  mutate(mod = mod/total,
         cyto = cyto/total) %>% 
  select(-total)

#### Find cluster vs isolated genes ####
map.all <- dat.all %>% 
  as.matrix() %>% 
  string_db$map(., "gene", removeUnmappedRows = TRUE) %>% 
  #collapse dups
  group_by(STRING_id) %>% 
  dplyr::summarise(gene = paste(unique(gene), collapse = " / "),
                   mod = mean(as.numeric(mod), na.rm=TRUE),
                   cyto = mean(as.numeric(cyto), na.rm=TRUE),
                   .groups="drop")
subgraph.all <- string_db$get_subnetwork(map.all$STRING_id)

isolated <- which(degree(subgraph.all)==0)
subgraph.all.clust <- delete.vertices(subgraph.all, isolated)

# Arrange metadata as in network
map.all.arrange <- map.all %>% 
  dplyr::filter(STRING_id %in% vertex_attr(subgraph.all.clust)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(subgraph.all.clust)$name)))

## Check order 
identical(vertex_attr(subgraph.all.clust)$name, map.all.arrange$STRING_id)
# Set attributes
##gene names
V(subgraph.all.clust)$symbol <- map.all.arrange$gene
##colors
vertex_attr(subgraph.all.clust)[["mod"]] <- unlist(map.all.arrange["mod"])
vertex_attr(subgraph.all.clust)[["cyto"]] <- unlist(map.all.arrange["cyto"])

#Layout
xy <- layout_with_fr(subgraph.all.clust) 

V(subgraph.all.clust)$x <- xy[, 1]
V(subgraph.all.clust)$y <- xy[, 2]

#plot
plot <- ggraph(subgraph.all.clust, layout= "manual", 
               x = V(subgraph.all.clust)$x, y = V(subgraph.all.clust)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey80") +
  scale_edge_width(range = c(0.2,1.5), name="STRING score") 

#Add nodes
plot + 
  geom_scatterpie(data=as_data_frame(subgraph.all.clust, "vertices"),
                  cols=sort(colnames(map.all.arrange)[-c(1:2)]), color=NA,
                  pie_scale = 0.7) +
  scale_fill_manual(values=color.vec, name="",
                    labels=c("cyto"="Cytokine protein",
                             "mod"="EOS02 module gene")) +
  geom_nodetext(aes(x = V(subgraph.all.clust)$x, y = V(subgraph.all.clust)$y,
                    label=V(subgraph.all.clust)$symbol), size=2) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")


#### Remove non-main cluster ####
#genes in small clusters
small.cluster <- c("KRT1","KRT6A","KRT81",
                   "CCT6A","VBP1",
                   "ADAMTS2", "ADAMTS14",
                   "THEM4","PYROXD2",
                   "MFAP4","SFTPD","SFTPA1","SFTPA2","SFTPC",
                   "JAZF1","CDC123",
                   "NAT8L","RIMKLA",
                   "ITGA7","COL4A1","COL9A3","COL19A1",
                   "PRRC2C","PIGC",
                   "SCT","GPR176",
                   "TNFRSF17","TNFRSF13B",
                   "PGD","TKTL1","G6PD","PGAM1","ENO1","PDHB","AK2","DLST",
                   "HSD17B6","UGT2B11",
                   "TMEM176A","TMEM176B",
                   "AFTPH","AP1S2","BLOC1S6",
                   "ANP32A","HMGB2",
                   "OTULIN","OTUD7B",
                   "DHDDS","PDSS1")

#Map genes to STRING
map <- dat.all %>% 
  filter(gene %notin% small.cluster) %>% 
  as.matrix() %>% 
  string_db$map(., "gene", removeUnmappedRows = TRUE) %>% 
  #collapse dups
  group_by(STRING_id) %>% 
  dplyr::summarise(gene = paste(unique(gene), collapse = " / "),
                   mod = mean(as.numeric(mod), na.rm=TRUE),
                   cyto = mean(as.numeric(cyto), na.rm=TRUE),
                   .groups="drop")
#Save
write_csv(map, file="publication/STRING.genes.csv")

#### Network creation ####
# Create igraph object 
subgraph <- string_db$get_subnetwork(map$STRING_id)

##Remove isolated nodes, keep large cluster
isolated <- which(degree(subgraph)==0)
subgraph.large.clust <- delete.vertices(subgraph, isolated)

#### Plot large cluster ####
current.graph <- subgraph.large.clust
current.map <- map

# Arrange metadata as in network
map.arrange <- current.map %>% 
  dplyr::filter(STRING_id %in% vertex_attr(current.graph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(current.graph)$name)))

## Check order 
identical(vertex_attr(current.graph)$name, map.arrange$STRING_id)
# Set attributes
##gene names
V(current.graph)$symbol <- map.arrange$gene
##colors
vertex_attr(current.graph)[["mod"]] <- unlist(map.arrange["mod"])
vertex_attr(current.graph)[["cyto"]] <- unlist(map.arrange["cyto"])

#Layout
set.seed(10)
xy <- layout_with_lgl(current.graph, cellsize=15) 

V(current.graph)$x <- xy[, 1]
V(current.graph)$y <- xy[, 2]

#plot
plot <- ggraph(current.graph, layout= "manual", 
               x = V(current.graph)$x, y = V(current.graph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey80") +
  scale_edge_width(range = c(0.1,1), name="STRING score") 

#Add nodes
plot.large.clust <- plot + 
    geom_scatterpie(data=as_data_frame(current.graph, "vertices"),
                    cols=sort(colnames(map.arrange)[-c(1:2)]), color=NA,
                    pie_scale = 0.5) +
    scale_fill_manual(values=color.vec, name="",
                      labels=c("cyto"="Cytokine protein",
                               "mod"="EOS02 module gene")) +
    geom_nodetext(aes(x = V(current.graph)$x, y = V(current.graph)$y,
                      label=V(current.graph)$symbol), size=2) +
    theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

#plot(plot.large.clust)
ggsave("publication/Fig2.STRING700_lrg.clust.pdf", plot.large.clust, 
       height=10, width=10)

#### Plot small clusters ####
current.graph <- subgraph.small.clust
current.map <- map.small.clust

# Arrange metadata as in network
map.arrange <- current.map %>% 
  dplyr::filter(STRING_id %in% vertex_attr(current.graph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(current.graph)$name)))

## Check order 
identical(vertex_attr(current.graph)$name, map.arrange$STRING_id)
# Set attributes
##gene names
V(current.graph)$symbol <- map.arrange$gene
##colors
vertex_attr(current.graph)[["mod"]] <- unlist(map.arrange["mod"])
vertex_attr(current.graph)[["cyto"]] <- unlist(map.arrange["cyto"])

#Layout
xy <- layout_with_fr(current.graph) 

V(current.graph)$x <- xy[, 1]
V(current.graph)$y <- xy[, 2]

#plot
plot <- ggraph(current.graph, layout= "manual", 
               x = V(current.graph)$x, y = V(current.graph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey80") +
  scale_edge_width(range = c(0.2,3), name="STRING score") 

#Add nodes
plot.small.clust <- plot + 
  geom_scatterpie(data=as_data_frame(current.graph, "vertices"),
                  cols=sort(colnames(map.arrange)[-c(1:2)]), color=NA,
                  pie_scale = 0.5) +
  scale_fill_manual(values=color.vec[2], name="",
                    labels=c("cyto"="Cytokine protein",
                             "mod"="EOS02 module gene")) +
  geom_nodetext(aes(x = V(current.graph)$x, y = V(current.graph)$y,
                    label=V(current.graph)$symbol), size=5) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

#plot(plot.small.clust)
ggsave("publication/Fig2.STRING400_sm.clust.pdf", plot.small.clust, height=30, width=30)

#### Plot isolated nodes ####
current.graph <- subgraph.isolate
current.map <- map

# Arrange metadata as in network
map.arrange <- current.map %>% 
  dplyr::filter(STRING_id %in% vertex_attr(current.graph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(current.graph)$name)))

## Check order 
identical(vertex_attr(current.graph)$name, map.arrange$STRING_id)
# Set attributes
##gene names
V(current.graph)$symbol <- map.arrange$gene
##colors
vertex_attr(current.graph)[["mod"]] <- unlist(map.arrange["mod"])
vertex_attr(current.graph)[["cyto"]] <- unlist(map.arrange["cyto"])

#Layout
#xy <- layout_with_fr(current.graph) 
xy <- layout_on_grid(current.graph) 

V(current.graph)$x <- xy[, 1]
V(current.graph)$y <- xy[, 2]

#plot
plot <- ggraph(current.graph, layout= "manual", 
               x = V(current.graph)$x, y = V(current.graph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey80") +
  scale_edge_width(range = c(0.2,3), name="STRING score") 

#Add nodes
plot.isolate <- plot + 
  geom_scatterpie(data=as_data_frame(current.graph, "vertices"),
                  cols=sort(colnames(map.arrange)[-c(1:2)]), color=NA,
                  pie_scale = 0.5) +
  scale_fill_manual(values=color.vec[2], name="",
                    labels=c("cyto"="Cytokine protein",
                             "mod"="EOS02 module gene")) +
  geom_nodetext(aes(x = V(current.graph)$x, y = V(current.graph)$y,
                    label=V(current.graph)$symbol), size=5) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

#plot(plot.isolate)
ggsave("publication/Fig2.STRING400_isolate.pdf", plot.isolate, height=30, width=30)


