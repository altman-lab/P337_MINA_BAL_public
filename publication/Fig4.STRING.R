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
                          score_threshold=400, input_directory="")

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

#### Find small cluster vs isolated genes ####
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

isolated <- which(degree(subgraph)==0)
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
                             "mod"="BAL EOS 02 gene")) +
  geom_nodetext(aes(x = V(subgraph.all.clust)$x, y = V(subgraph.all.clust)$y,
                    label=V(subgraph.all.clust)$symbol), size=2) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")


#### Map genes to STRING ####
#genes in small clusters
small.cluster <- c("KCNK17","KCNG2", "DHDDS","PDSS1", "PIGC","PRRC2C",
                   "PBX2","FAM222A", "MRO","MROH7", "VPREB3","C1ORF54",
                   "SDK2","KIAA1462", "FRMD6","FJX1", "HSPB9","SLC25A48","SLC25A53",
                   "PYROXD2","THEM4","SPIN4", "NAT8L","RIMKLA", "OTULIN","OTUD7B",
                   "VSIG1","CTSE")

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

map.small.clust <- dat.all %>% 
  filter(gene %in% small.cluster) %>% 
  as.matrix() %>% 
  string_db$map(., "gene", removeUnmappedRows = TRUE) %>% 
  #collapse dups
  group_by(STRING_id) %>% 
  dplyr::summarise(gene = paste(unique(gene), collapse = " / "),
                   mod = mean(as.numeric(mod), na.rm=TRUE),
                   cyto = mean(as.numeric(cyto), na.rm=TRUE),
                   .groups="drop")

#### Network creation ####
# Create igraph object 
subgraph <- string_db$get_subnetwork(map$STRING_id)
subgraph.small.clust <- string_db$get_subnetwork(map.small.clust$STRING_id)

##Remove isolated nodes, keep large cluster
isolated <- which(degree(subgraph)==0)
subgraph.large.clust <- delete.vertices(subgraph, isolated)

##remove clusters, keep isolated nodes
clusters <- which(degree(subgraph)>0)
subgraph.isolate <- delete.vertices(subgraph, clusters)

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
xy <- layout_with_lgl(current.graph) 

V(current.graph)$x <- xy[, 1]
V(current.graph)$y <- xy[, 2]

#plot
plot <- ggraph(current.graph, layout= "manual", 
               x = V(current.graph)$x, y = V(current.graph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey80") +
  scale_edge_width(range = c(0.2,1.5), name="STRING score") 

#Add nodes
plot.large.clust <- plot + 
    geom_scatterpie(data=as_data_frame(current.graph, "vertices"),
                    cols=sort(colnames(map.arrange)[-c(1:2)]), color=NA,
                    pie_scale = 0.5) +
    scale_fill_manual(values=color.vec, name="",
                      labels=c("cyto"="Cytokine protein",
                               "mod"="BAL EOS 02 gene")) +
    geom_nodetext(aes(x = V(current.graph)$x, y = V(current.graph)$y,
                      label=V(current.graph)$symbol), size=3) +
    theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

#plot(plot.large.clust)
ggsave("publication/Fig4.STRING400_lrg.clust.pdf", plot.large.clust, 
       height=14, width=12)

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
  scale_edge_width(range = c(0.2,1.5), name="STRING score") 

#Add nodes
plot.small.clust <- plot + 
  geom_scatterpie(data=as_data_frame(current.graph, "vertices"),
                  cols=sort(colnames(map.arrange)[-c(1:2)]), color=NA,
                  pie_scale = 0.7) +
  scale_fill_manual(values=color.vec[2], name="",
                    labels=c("cyto"="Cytokine protein",
                             "mod"="BAL EOS 02 gene")) +
  geom_nodetext(aes(x = V(current.graph)$x, y = V(current.graph)$y,
                    label=V(current.graph)$symbol), size=2) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

plot(plot.small.clust)
ggsave("publication/Fig4.STRING400_sm.clust.pdf", plot.small.clust, height=10, width=10)

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
  scale_edge_width(range = c(0.2,1.5), name="STRING score") 

#Add nodes
plot.isolate <- plot + 
  geom_scatterpie(data=as_data_frame(current.graph, "vertices"),
                  cols=sort(colnames(map.arrange)[-c(1:2)]), color=NA,
                  pie_scale = 0.7) +
  scale_fill_manual(values=color.vec[2], name="",
                    labels=c("cyto"="Cytokine protein",
                             "mod"="BAL EOS 02 gene")) +
  geom_nodetext(aes(x = V(current.graph)$x, y = V(current.graph)$y,
                    label=V(current.graph)$symbol), size=2) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

#plot(plot.isolate)
ggsave("publication/Fig4.STRING400_isolate.pdf", plot.isolate, height=10, width=10)


