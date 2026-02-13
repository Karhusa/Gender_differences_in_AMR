## 1.1 Load packages

```r
library(SummarizedExperiment)
library(vegan)
library(stats)
```
## 2.1 Upload TSE
```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

## 2.2 Inspect the TSE to create a new matrix from sample names and ARG genes (and abundances)
```r
counts <- assays(TSE)[[1]]  # usually the first assay
dim(counts)
head(counts[, 1:5])

genes <- rownames(TSE_object)  # often stored here
head(genes)

rowData(TSE_object)
```
Create the matrix:

* samples : rows
* genes : columns

```r


```

If your current matrix has rows = genes and columns = samples , transpose it:

```r
# counts: rows = genes, columns = samples
counts_matrix <- assays(TSE_object)[[1]]   # get the first assay
counts_matrix <- t(counts_matrix)         # transpose: now rows = samples, columns = genes
```
calculate the distance matrix

remember abundances

```r
# Convert counts to relative abundance per sample
counts_rel <- decostand(counts_matrix, method = "total") 


´´´


remember to sort out NA values (=zero counts)

```r
library(vegan)
# relative abundances:
dist_bc <- vegdis(arg_matrix, method = "bray")

as.matrix(dist_bc)[1:5, 1:5]
```




pcoa <- cmdscale(dist_bc, eig = TRUE, k = 2)

pcoa <- as.data.frame(reducedDim(tse, "PCoA"))



pcoa$age <- metada$age
PERMANOVA
adonis2
Dispersion
[12:51 PM]betadisper(dist_bc, meta$gender) (edited)
