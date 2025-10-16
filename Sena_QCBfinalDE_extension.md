QCB_455_final
================
Sena Cetin
2024-12-06

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (!require("locfit")) install.packages("locfit")
# BiocManager::install("edgeR", force=T)
library(BiocManager)
```

    ## Bioconductor version '3.18' is out-of-date; the current release version '3.21'
    ##   is available with R version '4.5'; see https://bioconductor.org/install

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
library(locfit)
```

    ## locfit 1.5-9.10   2024-06-24

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ purrr::none()   masks locfit::none()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(dplyr)
library(ShortRead)
```

    ## Loading required package: BiocGenerics
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## Loading required package: BiocParallel
    ## Loading required package: Biostrings
    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## Loading required package: XVector
    ## 
    ## Attaching package: 'XVector'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact
    ## 
    ## Loading required package: GenomeInfoDb
    ## 
    ## Attaching package: 'Biostrings'
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit
    ## 
    ## Loading required package: Rsamtools
    ## Loading required package: GenomicRanges
    ## Loading required package: GenomicAlignments
    ## Loading required package: SummarizedExperiment
    ## Loading required package: MatrixGenerics
    ## Loading required package: matrixStats
    ## 
    ## Attaching package: 'matrixStats'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count
    ## 
    ## 
    ## Attaching package: 'MatrixGenerics'
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars
    ## 
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'Biobase'
    ## 
    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians
    ## 
    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians
    ## 
    ## 
    ## Attaching package: 'GenomicAlignments'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     last
    ## 
    ## 
    ## Attaching package: 'ShortRead'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     id
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     compose
    ## 
    ## The following object is masked from 'package:tibble':
    ## 
    ##     view
    ## 
    ## The following objects are masked from 'package:locfit':
    ## 
    ##     left, right

``` r
library(ggplot2)
```

``` r
rna_seq = "/Users/senacetin/Downloads/GSE121279_rnaseq-counts-raw (1).txt"

TableOfCounts = read.table(rna_seq, header=TRUE)

foxp3_peaks <- read.csv(
  "/Users/senacetin/Documents/qcb455final/foxp3_peaks.csv"
)

# foxp1 peaks in treg
foxp1_peaks_treg <- read.csv(
  "/Users/senacetin/Documents/qcb455final/foxp1_peaks_Treg.csv"
)
# foxp1 peaks in tconv
foxp1_peaks_tconv <- read.csv(
  "/Users/senacetin/Documents/qcb455final/foxp1_peaks_Tconv.csv"
)
```

Heat map of differentially expressed (FDR-adjusted P\<.01) genes
associated with Foxp1 peaks clustered using k-means clustering.

## DE genes

``` r
library(readxl)
Samples <- read_excel("/Users/senacetin/Downloads/Samples.xlsx")

#read/use TableOfCounts
ALL_tpm <- TableOfCounts  

#define groups 
group <- factor(c(Samples$Group))

# create DGElist object
DEG <- DGEList(counts=ALL_tpm, group=group)

# filter genes w large counts for statistical analysis
keep <- filterByExpr(DEG)
DEG <- DEG[keep,,keep.lib.sizes=FALSE] # recomputes lib size

# nomrmalization
DEG <- calcNormFactors(DEG)

# group design matrix
design <- model.matrix(~0+group)
Groupnames <- unique(c(Samples$Treatment))
colnames(design) <- Groupnames

# fit for all data
DEG <- estimateDisp(DEG,design) 
fit <- glmQLFit(DEG,design)

# compare activated Treg samples
qlf.comp_a <- glmQLFTest(fit, contrast=c(-1,0,1,0))
top.comp_a <- topTags(qlf.comp_a) #Coefficient: -1*Foxp1_fl_aTreg 1*Foxp1_wt_aTreg
summary(dt<-decideTestsDGE(qlf.comp_a, adjust.method="BH", p.value = 0.01))
```

    ##        -1*Foxp1_fl_aTreg 1*Foxp1_wt_aTreg
    ## Down                                  844
    ## NotSig                              13711
    ## Up                                    801

``` r
#        -1*Foxp1_fl_aTreg 1*Foxp1_wt_aTreg
# Down                                  801
# NotSig                              13711
# Up                                    844

#export results
write.csv(qlf.comp_a, "Foxp1_fl_aTreg_vs_Foxp1_wt_aTreg.csv") 

#############

# compare naive Treg samples
qlf.comp_n<- glmQLFTest(fit, contrast=c(0,-1,0,1))
topTags(qlf.comp_n) #Coefficient: Coefficient: -1*Foxp1_fl_nTreg 1*Foxp1_wt_nTreg 
```

    ## Coefficient:  -1*Foxp1_fl_nTreg 1*Foxp1_wt_nTreg 
    ##             logFC   logCPM         F       PValue          FDR
    ## Qpct    5.6897906 5.120420 1139.6532 3.597795e-14 5.524774e-10
    ## Ptprf   3.8588754 5.861713  501.0082 7.340518e-12 5.636050e-08
    ## Slamf6  1.2268833 9.150537  430.0310 1.954171e-11 1.000275e-07
    ## Tkt     1.2236868 8.395753  360.9865 5.969142e-11 2.291554e-07
    ## Cd1d1   1.5795513 7.102555  341.5503 8.486160e-11 2.606269e-07
    ## Slc9a7 -2.8118890 4.470594  284.5475 2.695652e-10 6.899071e-07
    ## Prkcd   0.9533412 7.213723  251.1748 5.908473e-10 1.259268e-06
    ## Kcnj16  5.2967535 3.965943  242.9822 7.273553e-10 1.259268e-06
    ## Bicd1  -3.5210018 3.734504  242.4169 7.380443e-10 1.259268e-06
    ## Endod1 -1.8589863 6.599338  226.6752 1.123028e-09 1.724522e-06

``` r
summary(dt<-decideTestsDGE(qlf.comp_n, adjust.method="BH", p.value = 0.01))
```

    ##        -1*Foxp1_fl_nTreg 1*Foxp1_wt_nTreg
    ## Down                                  351
    ## NotSig                              14766
    ## Up                                    239

``` r
#        -1*Foxp1_fl_nTreg 1*Foxp1_wt_nTreg
# Down                                  239
# NotSig                              14766
# Up                                    351

# export results
write.csv(qlf.comp_n, "Foxp1_fl_nTreg_vs_Foxp1_wt_nTreg.csv") 

###################
```

Heatmaps

``` r
#install.packages('pheatmap')
library(pheatmap)

counts.in <- TableOfCounts

counts.metadata <- data.frame(
  dataset= c(colnames(counts.in)), 
  Treatment = c(Samples$Treatment),
  stringsAsFactors = FALSE)

group_HM <- counts.metadata$Treatment # sample #

# create DGEList object 
y_HM <- DGEList(counts=counts.in, genes=row.names.data.frame(counts.in), group=group_HM)

# Low count removal -- We don't know how they filtered for low counts 
keep_HM <- rowSums(cpm(y_HM)>1) >= 1
table(keep_HM)
```

    ## keep_HM
    ## FALSE  TRUE 
    ## 32581 14301

``` r
y_HM <- y_HM[keep_HM, , keep.lib.sizes=FALSE]

# again, don't know what method used to normalize lib size
y_HM <- calcNormFactors(y_HM, method="TMM") 

# new table of counts
counts_HM <- as.matrix(y_HM$counts)

#normalization using cpm here
logCPM_HM <- cpm(counts_HM, prior.count=1, log=TRUE) 

# z-score matrix for heatmap
ZScore_HM <- t(scale(t(logCPM_HM)))
ZScore_HM <- as.data.frame(ZScore_HM)

# import DEGs for naive and activated samples 
comp_n <- read.csv(
  "~/Downloads/Foxp1_fl_nTreg_vs_Foxp1_wt_nTreg.csv", 
  col.names = c("Symbol", "logFC", "logCPM", "F", "PValue")
  )
comp_a <- read.csv("~/Downloads/Foxp1_fl_aTreg_vs_Foxp1_wt_aTreg.csv", 
  col.names = c("Symbol", "logFC", "logCPM", "F", "PValue")
  )

Sig_all_naive <- subset.data.frame(comp_n, PValue<0.01) #All significant naive DEG
Sig_all_activated <- subset.data.frame(comp_a, PValue<0.01) #All significant activated DEG
Sig_all <- rbind(Sig_all_naive, Sig_all_activated) #All sig DEGs for both comparisons

# filter for foxp1 bound
Sig_all_f1 <- Sig_all %>% 
  filter(Symbol %in% foxp1_peaks_treg$gene) # 1015 genes 

# filter for unique genes 
Sig_all_unique <- Sig_all_f1[!duplicated(Sig_all_f1$Symbol), ] # 823 unique genes

# pull list of gene symbols
allSymbols <- c(Sig_all_unique$Symbol)

# get zscores and store as matrix
sig.zscore <- ZScore_HM[allSymbols,] 
sig.zscore.mat <- as.matrix(sig.zscore)
sig.zscore.mat <- sig.zscore.mat[complete.cases(sig.zscore.mat),] # no missing vals

# pull annotations from metadata for heatmap
heat.annotation <- data.frame(counts.metadata[,2])
colnames(heat.annotation) <- "Treatment"
row.names(heat.annotation) <- counts.metadata[,1]
nrow(sig.zscore.mat)
```

    ## [1] 823

``` r
colnames(sig.zscore.mat) <- c(counts.metadata[,1])

# assign annotation colors from Samples
ann.colors <- list(Treatment= c(
  `Foxp1_fl_aTreg` = Samples$Group_color[3],
  `Foxp1_fl_nTreg` = Samples$Group_color[6], 
  `Foxp1_wt_aTreg` = Samples$Group_color[9],
  `Foxp1_wt_nTreg` = Samples$Group_color[12]))

# create table for annotations for each gene if they are foxp1 bound, foxp3 bound,
# up regulated and down regulated for foxp1 deficient activated and naive tregs
generow_annotations <- data.frame(Symbol = allSymbols) %>% 
  mutate(Foxp1 = as.factor(if_else(Symbol %in% foxp1_peaks_treg$gene, 1, NA))) %>% 
  mutate(Foxp3 = as.factor(if_else(Symbol %in% foxp3_peaks$gene, 1, NA))) %>% 
  # find up and down-regulated in foxp1 deficient NAIVE treg cells
  left_join(Sig_all_naive %>%  # Add logFC_naive
              dplyr::select(Symbol, logFC) %>% 
              dplyr::rename(logFC_naive = logFC),
              by = "Symbol") %>%
  # annotated as up or down regulated based on logFC
  mutate(expression_naive = case_when(
    logFC_naive >= 0 ~ "Upregulated", 
    logFC_naive < 0 ~ "Downregulated",
    TRUE ~ NA_character_
  )) %>% 
  # find up and down-regulated in foxp1 deficient ACTIVATED treg cells
  left_join(Sig_all_activated %>%  # Add logFC_naive
            dplyr::select(Symbol, logFC) %>% 
            dplyr::rename(logFC_activated = logFC),
            by = "Symbol") %>%
  # annotated as up or down regulated based on logFC
  mutate(expression_activated = case_when(
    logFC_activated >= 0 ~ "Upregulated", 
    logFC_activated < 0 ~ "Downregulated",
    TRUE ~ NA_character_
  )) %>%
  select(-logFC_naive) %>% # remove logFC cols 
  select(-logFC_activated)

row_annotations <- data.frame(
  Symbol = generow_annotations$Symbol,
  "Foxp1 Bound" = generow_annotations$Foxp1,
  "Foxp3 Bound" = generow_annotations$Foxp3,
  "Foxp1 deficient nTreg" = generow_annotations$expression_naive,
  "Foxp1 deficient aTreg" = generow_annotations$expression_activated
) %>% 
  column_to_rownames(var = "Symbol")  # Convert `Symbol` to row names

# Check alignment of row names
stopifnot(all(rownames(sig.zscore.mat) %in% rownames(row_annotations)))

dist_rows <- dist(sig.zscore.mat)
tree_row  <- hclust(dist_rows, method = "complete")
row_clusters <- cutree(tree_row, 6)

heatmap_all_hierarch <- pheatmap(sig.zscore.mat,
                        annotation_col = heat.annotation, 
                        cluster_cols = TRUE,
                        cluster_rows = tree_row,
                        clustering_method = "complete", # comp
                        main = "Hierarchical Clustered DEGs Associated With Foxp1 Peaks",
                        annotation_colors = ann.colors,
                        show_colnames = TRUE,
                        show_rownames = FALSE,
                        annotation_row = row_annotations,
                        #labels_row = as.expression(newnames),
                        #fontsize_row = 6,
                        cutree_rows = 6,
                        gaps_row = cumsum(table(row_clusters))
                        )
```

<img src="Sena_QCBfinalDE_extension_files/figure-gfm/counts.in for heatmap-1.png" width="100%" style="display: block; margin: auto;" />

``` r
#plot_width <- 3
#plot_height <- 4

# Save the plot using ggsave
# ggsave("/Users/senacetin/Downloads/heatmap_all_hierarch.png",
#       plot = heatmap_all_hierarch$gtable,
#       width = plot_width, height = plot_height, limitsize = FALSE)

# mutate column of cluster number to all significant foxp1 bound DEGs
clusters_hier <- Sig_all_f1 %>% 
  mutate(Cluster = row_clusters[Symbol])

# Export DEG list with associated cluster as CSV for downstream analysis
write.csv(clusters_hier, "clusters_hier.csv", row.names = TRUE)
```

``` r
#install.packages("dyndimred")
library(dyndimred)
#install.packages("Rtsne")
library(Rtsne)
# install.packages('umap')
library(umap)

et <- exactTest(DEG)
topTags(et)
```

    ## Comparison of groups:  2-1 
    ##              logFC   logCPM        PValue           FDR
    ## Asb2     -8.357145 4.594323 2.141940e-257 3.289163e-253
    ## Fgl2     -5.127660 5.155141 8.969740e-154 6.886966e-150
    ## Tnfrsf1b -1.905273 8.806154 9.007792e-136 4.610789e-132
    ## Itgb1    -3.048872 9.769597 1.652408e-128 6.343595e-125
    ## Ttc39c   -6.646469 4.838947 2.266932e-119 6.962201e-116
    ## Il9r     -5.214260 4.262875 1.806424e-112 4.623242e-109
    ## Cd83     -2.274412 7.451991 1.977568e-101  4.338219e-98
    ## Fam129a  -2.178100 7.109958 6.352732e-100  1.219407e-96
    ## Gpr68    -3.331760 6.295899  7.590647e-98  1.295133e-94
    ## Abhd4    -2.944287 4.507520  4.058673e-90  6.232498e-87

``` r
DEG <- estimateDisp(DEG,design)
DEG <- estimateGLMTagwiseDisp(DEG, design)
logcpm <- cpm(DEG, log=TRUE)

set.seed(123) 
t_logcpm <- t(logcpm)
all_labels=c(Samples$Name)


####
DEG <- estimateDisp(DEG, design, robust=TRUE) 
DEG$common.dispersion
```

    ## [1] 0.02660844

``` r
sqrt(DEG$common.dispersion)
```

    ## [1] 0.1631209

``` r
plotBCV(DEG)
```

![](Sena_QCBfinalDE_extension_files/figure-gfm/pca_tsne_Umap_all-1.png)<!-- -->

``` r
######

GroupName <- unique(Samples$Treatment) 

#PCA 
pca <- prcomp(t_logcpm, scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], 
     pch = 19, cex=2, 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA of RNA-seq Samples",
     col=Samples$Group_color, asp=0)
legend("topright", 
       legend = GroupName, 
       col = unique(Samples$Group_color), 
       pch = 19, cex = 1, xpd = TRUE, box.lwd = 0, 
       bg = "transparent")
```

![](Sena_QCBfinalDE_extension_files/figure-gfm/pca_tsne_Umap_all-2.png)<!-- -->

``` r
summary(pca)
```

    ## Importance of components:
    ##                            PC1     PC2      PC3      PC4      PC5     PC6
    ## Standard deviation     70.0423 57.2709 35.47830 32.60229 31.50927 28.6900
    ## Proportion of Variance  0.3195  0.2136  0.08197  0.06922  0.06465  0.0536
    ## Cumulative Proportion   0.3195  0.5331  0.61504  0.68426  0.74891  0.8025
    ##                             PC7     PC8     PC9     PC10     PC11      PC12
    ## Standard deviation     27.02702 25.6055 24.7231 23.39482 22.08847 4.708e-13
    ## Proportion of Variance  0.04757  0.0427  0.0398  0.03564  0.03177 0.000e+00
    ## Cumulative Proportion   0.85009  0.8928  0.9326  0.96823  1.00000 1.000e+00

``` r
#tSNE 
tsne_out <- Rtsne(t_logcpm, pca=TRUE, perplexity=3, theta=0.0, check_duplicates=TRUE)

GroupName <- unique(Samples$Treatment) 

# Create the scatter plot with legend
plot(tsne_out$Y, 
     pch = 19, cex = 2, 
     xlab = "tSNE1", ylab = "tSNE2",
     main = "tSNE of RNA-seq Samples",
     col = Samples$Group_color, asp = 0)
legend("topleft", 
       legend = GroupName, 
       col = unique(Samples$Group_color), 
       pch = 19, cex = 1, xpd = TRUE, box.lwd = 0, 
       bg = "transparent")
```

![](Sena_QCBfinalDE_extension_files/figure-gfm/pca_tsne_Umap_all-3.png)<!-- -->

``` r
# UMAP
umap_out <- umap(t_logcpm,
                method = "naive",
                n_neighbors = 2, # each group has 1+2 neighbors
                n_components = 4, # four groups
                metric = "euclidean",
                preserve.seed = TRUE,
                )
```

    ## Warning: failed creating initial embedding; using random embedding instead
    ## Warning: failed creating initial embedding; using random embedding instead
    ## Warning: failed creating initial embedding; using random embedding instead
    ## Warning: failed creating initial embedding; using random embedding instead

``` r
plot(umap_out$layout, 
     pch = 19, cex = 2, 
     xlab = "UMAP1", ylab = "UMAP2", 
     main = "UMAP of RNA-seq Samples",
     col = Samples$Group_color, asp = 0)
legend("topright", 
       legend = GroupName, 
       col = unique(Samples$Group_color), 
       pch = 19, cex = 1, xpd = TRUE, box.lwd = 0, bg = "transparent")
```

![](Sena_QCBfinalDE_extension_files/figure-gfm/pca_tsne_Umap_all-4.png)<!-- -->
