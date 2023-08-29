library(tidyverse)
library(patchwork)
library(readxl)
library(ggrepel)

#### RNAseq data ####
load("data_clean/P337_BAL_data.RData")
load("data_clean/P337_BAL_module_data.RData")

#### Gene PCA ####
PCA <- as.data.frame(dat.BAL.abund.norm.voom$E) %>% 
  t() %>% 
  prcomp()

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1], digits=3)*100, "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2], digits=3)*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.BAL.abund.norm.voom$targets, by="libID")

PCA_gene <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit, shape=visit), size=2) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, title="Genes") +
  coord_fixed(ratio=1) +
  theme(legend.position = "none",
        plot.title = element_text(size=10)) +
  scale_color_manual(values = c(V4 = "#74add1", V5 = "#FF8679")) 

#### Module PCA ####
PCA <- as.data.frame(mod.voom) %>% 
  column_to_rownames("module") %>% 
  t() %>% 
  prcomp()

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1], digits=3)*100, "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2], digits=3)*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.BAL.abund.norm.voom$targets, by="libID")

PCA_mod <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit, shape=visit), size=2) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, title="Modules") +
  coord_fixed(ratio=1)+
  scale_colour_manual(values = c(V4 = "#74add1", V5 = "#FF8679"),
                      name = "SBP-Ag",
                      labels = c("Pre","Post")) +   
  scale_shape_manual(name = "SBP-Ag",
                     labels = c("Pre","Post"),
                     values = c(19,17)) +
  theme(plot.title = element_text(size=10))

#### Cytokine PCA ####
plex <- read_csv("data_raw/addtl.data/P337_BAL.multiplex.csv") %>% 
  rename(donorID=ptID, FGF2_V4=TGF2_V4) %>% 
  pivot_longer(-donorID) %>% 
  #Log10 transform
  mutate(value = log10(value)) %>% 
  separate(name, into=c("rowname","visit"), sep="_") %>% 
  mutate(name =paste(donorID, visit, sep="_")) %>% 
  select(-donorID, -visit) %>% 
  pivot_wider() %>% 
  column_to_rownames() %>% 
  t()

PCA <- plex %>% 
  prcomp()

PC1.label <- paste("PC1 (", round(summary(PCA)$importance[2,1], digits=3)*100, "%)", sep="")
PC2.label <-paste("PC2 (", round(summary(PCA)$importance[2,2], digits=3)*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("sampID") %>%
  # Select PCs
  dplyr::select(sampID, PC1:PC3) %>% 
  # Add  metadata
  separate(sampID, into=c("donorID","visit"), sep="_", remove = FALSE)

PCA_prot <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_line(aes(group=donorID), color="grey") +
  geom_point(aes(color=visit, shape=visit), size=2) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, title="Proteins") +
  coord_fixed(ratio=1) +
  theme(legend.position = "none",
        plot.title = element_text(size=10))+
  scale_color_manual(values = c(V4 = "#74add1", V5 = "#FF8679")) 

#### Gene volcano ####
gene_visit <- read_csv("results/gene_level/P337_BAL_gene_visit.csv") %>% 
  rename(estimate=logFC, FDR=adj.P.Val) %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  rename(gene=hgnc_symbol) %>% 
  mutate(variable = recode(group, "visit"="Post - Pre SBP-Ag"))

gene_EOS <- read_csv("results/gene_level/P337_BAL_gene_EOS.csv")%>% 
  rename(estimate=logFC, FDR=adj.P.Val) %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  rename(gene=hgnc_symbol) %>% 
  mutate(variable = recode(group, "EOS.pct"="Eosinophil cell %"))

gene_PMN <- read_csv("results/gene_level/P337_BAL_gene_NEUT.csv")%>% 
  rename(estimate=logFC, FDR=adj.P.Val) %>% 
  left_join(dat.BAL.abund.norm.voom$genes) %>% 
  rename(gene=hgnc_symbol) %>% 
  mutate(variable = recode(group, "NEUT.pct"="Neutrophil cell %"))

model_result <- bind_rows(gene_visit, gene_EOS, gene_PMN) %>% 
  filter(variable != "(Intercept)") %>% 
  mutate(variable = fct_relevel(variable, "Post - Pre SBP-Ag", after = 0L))

top3up <- model_result %>%
  filter(FDR < 0.3 & estimate > 0) %>% 
  group_by(variable) %>% 
  slice_min(FDR, n=3) %>% 
  mutate(lab = gene)

top3dn <- model_result %>%
  filter(FDR < 0.3 & estimate < 0) %>% 
  group_by(variable) %>% 
  slice_min(FDR, n=3) %>% 
  mutate(lab = gene)

model_format <- model_result %>% 
  mutate(col.group = case_when(FDR < 0.3 & estimate < 0 ~ "down", 
                               FDR < 0.3 & estimate > 0 ~ "up", 
                               TRUE ~ "NS")) %>%
  
  full_join(bind_rows(top3up,top3dn)) %>% 
  mutate(col.group = factor(col.group, levels = c("down", "up", "NS"))) %>% 
  arrange(desc(col.group))

volc_gene <- ggplot(model_format, aes(x = estimate, y = -log10(FDR))) + 
  geom_point(aes(color = col.group)) +
  theme_minimal() + 
  labs(y = "-log10( FDR )", x = "Log2 fold change", color = "FDR < 0.3") + 
  facet_wrap(~variable, scales = "free_x") + 
  scale_color_manual(values = c(down = "#74add1", NS = "grey", up = "#FF8679"), 
                     na.value = "grey") +
  # geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_hline(yintercept = -log10(0.3), lty = "dashed") +
  geom_text_repel(aes(label = lab), direction = "both", 
                  min.segment.length = ggplot2::unit(0, "lines"), 
                  show.legend = FALSE, max.overlaps = Inf,
                  size=3) +
  theme(panel.border = element_rect(fill=NA))

#### Protein volcano ####
prot_visit <- read_csv("results/protein_level/P337_BAL_protein_visit.csv")%>% 
  rename(estimate=logFC, FDR=adj.P.Val) %>% 
  rename(gene=geneName) %>% 
  mutate(variable = recode(group, "visit"="Post - Pre SBP-Ag")) %>% 
  filter(variable != "(Intercept)")

top3up <- prot_visit %>%
  filter(FDR < 0.3 & estimate > 0) %>% 
  group_by(variable) %>% 
  slice_min(FDR, n=3) %>% 
  mutate(lab = gene)

top3dn <- prot_visit %>%
  filter(FDR < 0.3 & estimate < 0) %>% 
  group_by(variable) %>% 
  slice_min(FDR, n=3) %>% 
  mutate(lab = gene)

prot_visit_format <- prot_visit %>% 
  mutate(col.group = case_when(FDR < 0.3 & estimate < 0 ~ "down", 
                               FDR < 0.3 & estimate > 0 ~ "up", 
                               TRUE ~ "NS")) %>%
  
  full_join(bind_rows(top3up,top3dn)) %>% 
  mutate(col.group = factor(col.group, levels = c("down", "up", "NS"))) %>% 
  arrange(desc(col.group))

volc_prot <- ggplot(prot_visit_format, aes(x = estimate, y = -log10(FDR))) + 
  geom_point(aes(color = col.group)) +
  theme_minimal() + 
  labs(y = "-log10( FDR )", x = "Log2 fold change", color = "FDR < 0.3") + 
  facet_wrap(~variable, scales = "free_x") + 
  scale_color_manual(values = c(down = "#74add1", NS = "grey", up = "#FF8679"), 
                     na.value = "grey") +
  # geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_hline(yintercept = -log10(0.3), lty = "dashed") +
  geom_text_repel(aes(label = lab), direction = "both", 
                  min.segment.length = ggplot2::unit(0, "lines"), 
                  show.legend = FALSE, max.overlaps = Inf,
                  size=3) +
  theme(panel.border = element_rect(fill=NA), legend.position = "none")

#### Save ####
lo <- "
ACE
BDD
"
plot_all <- PCA_prot + volc_prot + PCA_gene + volc_gene + PCA_mod +
  plot_layout(design = lo) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
# plot_all

ggsave(plot_all, filename="publication/Fig1.PCA.volcano.png", height=6, width=10)
ggsave(plot_all, filename="publication/Fig1.PCA.volcano.pdf", height=6, width=10)

