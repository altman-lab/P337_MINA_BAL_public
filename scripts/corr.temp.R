## Pearson correlation

```{r}
corr.all <- full_join(counts.mod.delta, neuro.delta, by = "donorID") %>% 
  full_join(asthma.delta, by = "donorID") %>% 
  full_join(plex.delta, by="donorID")
```

For variables that are retained in (S)PLS, run Pearson correlations

```{r warning=FALSE, message=FALSE}
#Subset data to PLS variables
splsX.var <- as.data.frame(result.spls$loadings$X) %>% 
  filter(comp1 != 0 | comp2 != 0) %>% 
  rownames_to_column() %>% 
  distinct(rowname) %>% unlist(use.names = FALSE)

splsY.var <- as.data.frame(result.spls$loadings$Y) %>% 
  filter(comp1 != 0 | comp2 != 0) %>% 
  rownames_to_column() %>% 
  distinct(rowname) %>% unlist(use.names = FALSE)

cor.sub <- corr.all %>% 
  dplyr::select(donorID, all_of(splsX.var), all_of(splsY.var))

corr.fxn(counts=cor.sub, 
         metadata=cor.sub, 
         match.by="donorID", 
         corr="pearson", 
         basename="P337_delta", 
         heatmap=TRUE, dotplot=FALSE)
```

Load results.

```{r message=FALSE}
#R^2 values
R <- read_csv("results/correlation/R_pearson_P337_delta.csv")
#P-values
P <- read_csv("results/correlation/P_pearson_P337_delta.csv")
```

#### Heatmaps

All correlations

```{r echo=FALSE}
R.format <- R %>% 
  mutate(variable = gsub(".multiplex|.pct", "", variable)) %>% 
  rename_all(~gsub(".multiplex|.pct", "", .)) %>% 
  column_to_rownames("variable") %>% 
  as.matrix()
#diagonal to NA
diag(R.format) <- NA
```

```{r echo=FALSE}
# Calculate FDR
P.lower <- P %>% 
  column_to_rownames("variable")
P.lower[upper.tri(P.lower, diag=FALSE)] <- NA

FDR.lower <-as.data.frame(P.lower) %>%
  rownames_to_column("variable") %>% 
  pivot_longer(-variable) %>% 
  #FDR groups
  mutate(value = p.adjust(value, method="BH")) 

P.upper <- P %>% 
  column_to_rownames("variable")
P.upper[lower.tri(P.upper, diag=FALSE)] <- NA

FDR.upper <-  P.upper %>%
  rownames_to_column("variable") %>% 
  pivot_longer(-variable) %>% 
  mutate(value = p.adjust(value, method="BH"))

FDR.format <- bind_rows(FDR.lower, FDR.upper) %>% 
  distinct() %>% 
  drop_na(value) %>% 
  arrange(match(name, P$variable)) %>% 
  pivot_wider() %>% 
  arrange(match(variable, P$variable)) %>% 
  mutate(variable = gsub(".multiplex|.pct", "", variable)) %>% 
  rename_all(~gsub(".multiplex|.pct", "", .)) %>% 
  column_to_rownames("variable") %>% 
  as.matrix()
```

```{r include=FALSE}
#Check order
identical(rownames(R.format), rownames(FDR.format))
identical(colnames(R.format), colnames(FDR.format))
```

```{r echo=FALSE, message=FALSE, fig.width=8.5, fig.height=9}
corrplot(R.format,
         method="color", order="original", type="upper", diag=TRUE,
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90,
         #Change color legend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = FDR.format,
         insig = "label_sig",
         sig.level = c(0.01, 0.05), 
         pch.cex=0.9, pch.col="white",
         na.label = " ",
         #Add title
         title="FDR *<0.05 **<0.01", 
         mar = c(0,0,2,0))

#### Save ####

pdf(file = "figs/correlation/P337_delta/all_pearson_P337_delta_format.pdf", 
    height=10, width=10)

corrplot(R.format,
         method="color", order="original", type="upper", diag=TRUE,
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90,
         #Change color legend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = FDR.format,
         insig = "label_sig",
         sig.level = c(0.01, 0.05), 
         pch.cex=0.9, pch.col="white",
         na.label = " ",
         #Add title
         title="FDR *<0.05 **<0.01", 
         mar = c(0,0,2,0))

dev.off()
```

Significant correlations to fMRI

```{r echo=FALSE}
#Significant for at least 1 fMRI
to.keepX <- as.data.frame(FDR.format) %>% 
  rownames_to_column("variable") %>% 
  dplyr::select(variable, R_vA_insula:LPR_ACC_EOS) %>% 
  filter(grepl("BAL_", variable) | !grepl("insula|amgy|amyg|ACC", variable)) %>% 
  filter_at(vars(-variable), any_vars(.<0.05))
#Significant for at least 1 mod or cytokine
to.keepY <- as.data.frame(FDR.format) %>% 
  rownames_to_column("variable") %>% 
  dplyr::select(-c(R_vA_insula:LPR_ACC_EOS)) %>% 
  filter(grepl("insula|amgy|amyg|ACC", variable)) %>% 
  filter_at(vars(-variable), any_vars(.<0.05))

R.signif <- as.data.frame(R.format) %>% 
  rownames_to_column("variable") %>% 
  filter(variable %in% to.keepX$variable | 
           variable %in% to.keepY$variable) %>% 
  dplyr::select(variable, all_of(to.keepX$variable), 
                all_of(to.keepY$variable)) %>% 
  column_to_rownames("variable") %>% 
  as.matrix()

FDR.signif <- as.data.frame(FDR.format) %>% 
  rownames_to_column("variable") %>% 
  filter(variable %in% to.keepX$variable | 
           variable %in% to.keepY$variable) %>% 
  dplyr::select(variable, all_of(to.keepX$variable), 
                all_of(to.keepY$variable)) %>% 
  column_to_rownames("variable") %>% 
  as.matrix()
```

```{r include=FALSE}
#Check order
identical(rownames(R.signif), rownames(FDR.signif))
identical(colnames(R.signif), colnames(FDR.signif))
```

```{r echo=FALSE, message=FALSE, fig.width=8.5, fig.height=15}
corrplot(R.signif,
         method="color", order="original", type="upper", diag=TRUE,
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90,
         #Change color legend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = FDR.signif,
         insig = "label_sig",
         sig.level = c(0.01, 0.05), 
         pch.cex=0.9, pch.col="white",
         na.label = " ",
         #Add title
         title="FDR *<0.05 **<0.01", 
         mar = c(0,0,2,0))

#### Save ####

pdf(file = "figs/correlation/P337_delta/signif_pearson_P337_delta_format.pdf", 
    height=5, width=5)

corrplot(R.signif,
         method="color", order="original", type="upper", diag=TRUE,
         #Reverse color order from default
         col=colorRampPalette(c("darkblue", "white",
                                "darkred"))(20),
         #Change labels
         tl.col="black", tl.srt=90,
         #Change color legend
         cl.pos="b", cl.length = 5,
         #Add significant labels
         p.mat = FDR.signif,
         insig = "label_sig",
         sig.level = c(0.01, 0.05), 
         pch.cex=0.9, pch.col="white",
         na.label = " ",
         #Add title
         title="FDR *<0.05 **<0.01", 
         mar = c(0,0,2,0))

dev.off()
```

