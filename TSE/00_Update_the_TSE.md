Update the original TSE.rds with new coldata table

## 1. Load the TSE object

```r

library(SummarizedExperiment)
library(S4Vectors)

# Load TSE
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)

# convert age_category to character to avoid factor warnings
if ("age_category" %in% colnames(colData(TSE))) {
  colData(TSE)$age_category <- as.character(colData(TSE)$age_category)
}

# Load new colData TSV
coldata_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/colData_TSE.tsv"
colData_new <- read.delim(coldata_path, header = TRUE, stringsAsFactors = FALSE)

# Extract sample IDs
acc_tse <- colData(TSE)$acc
acc_tsv <- colData_new$acc

# Find matching samples
common_acc <- intersect(acc_tse, acc_tsv)
length(common_acc)  # How many will be updated

# For each column in the TSV (except 'acc'), update TSE in-place
cols_to_update <- setdiff(colnames(colData_new), "acc")

for (colname in cols_to_update) {
  # If column doesn't exist in TSE yet, create it
  if (!colname %in% colnames(colData(TSE))) {
    TSE[[colname]] <- NA
  }
  # Update only matching rows
  TSE[[colname]][match(common_acc, acc_tse)] <- colData_new[[colname]][match(common_acc, acc_tsv)]
}

# Save updated TSE
saveRDS(TSE, tse_path)

```
## 2. Save the new colData as a TSV file
```
coldata_df <- as.data.frame(colData(TSE))

write.table(
  coldata_df,
  file = "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/colData_TSE_updated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```


