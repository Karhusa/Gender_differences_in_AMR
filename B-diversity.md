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

rowData(TSE_object)
```
## 2.3 Create the matrix:

* samples : rows
* genes : columns

```r
# 1. Get counts
counts <- assays(TSE)[[1]]   # 1406 genes x 60996 columns

# 2. Get sample names (accessions)
sample_ids <- colData(TSE)$acc

# 3. Create new matrix: rows = samples, columns = genes
counts_matrix <- t(counts)          # transpose: now 60996 x 1406
rownames(counts_matrix) <- sample_ids
colnames(counts_matrix) <- rownames(TSE)

# 4. Sanity check
dim(counts_matrix)
head(counts_matrix[, 1:5])
head(rownames(counts_matrix))
```

## 2.4 Save the matrix

```
saveRDS(counts_matrix, file = "counts_matrix.rds")

#counts_matrix <- readRDS("counts_matrix.rds")

```
## 2.5 Check if values are as relative abundances
```r
row_sums <- rowSums(counts_matrix)

summary(row_sums)
head(row_sums)


```
# summary(row_sums)
#     Min.   1st Qu.    Median Mean   3rd Qu.      Max. 
#    50.02    316.68    509.45 713.97    827.73 115515.58 
     
# head(row_sums)
# DRR070899 DRR070900 DRR070901 DRR070902 DRR070903 DRR070904 
#  760.1088  581.5741  600.2330  647.9989  567.5491  767.9853 


## 2.6 Convert counts to relative abundance per sample

# Reserve at least 16 GB

```r
counts_rel <- decostand(counts_matrix, method = "total") 


´´´

calculate the distance matrix
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
