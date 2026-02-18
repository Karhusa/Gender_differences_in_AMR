# Workflow: Matching TSE Samples to Metalog / SRA Data

Input files:
* TSE-object (TSE-filtered) from study (Mahkameh et. al, 2025 [https://www.nature.com/articles/s41522-025-00715-9])
* Metalog Samples (human_sample_list.csv) from (Metalog [https://metalog.embl.de/explore/human])
* SRA samples (SRA_metadata_with_biosample.txt) from Katariina Pärnänen

---

## 1. Load and Inspect TSE Object

```r
library(tibble)

TSE_filtered <- readRDS("~/Downloads/TSE_filtered.rds")

df_col <- as_tibble(colData(TSE_filtered))
View(df_col)
```
Goal: Check which accession numbers are present in the TSE object.

---

## 2. Download Sample List from Metalog and Inspect Accession Numbers

File: human_sample_list.csv

Columns:
1. Study code
2. ena_ers_sample_id (our main interest)
3. Sample alias

Use Unix tools for large files:

```bash
cut -d',' -f1 human_sample_list.csv | sort -u > studies.txt

tail -n +2 human_sample_list.csv | cut -d',' -f2 | sed 's/[0-9]*//g' | sort -u > unique_sample_names.txt
```
**Example prefixes from unique_sample_names:** ERR, ERS, SAMD, SAMEA, SAMN, SRR, SRS, SRX

---

## 3. Create Separate Tables for Each Prefix

```bash
awk -F',' 'NR==1 || $2 ~ /^ERR/ {print}' human_sample_list.csv > err_ids.csv
awk -F',' 'NR==1 || $2 ~ /^ERS/ {print}' human_sample_list.csv > ers_ids.csv
awk -F',' 'NR==1 || $2 ~ /^SAMD/ {print}' human_sample_list.csv > samd_ids.csv
awk -F',' 'NR==1 || $2 ~ /^SAMEA/ {print}' human_sample_list.csv > samea_ids.csv
awk -F',' 'NR==1 || $2 ~ /^SAMN/ {p
```
Goal: Smaller files for easier processing in R.

---

## 4. Load Split Sample Lists into R
```r
library(readr)

err_ids <- read_csv("Gradu_AMR/err_ids.csv")
ers_ids <- read_csv("Gradu_AMR/ers_ids.csv")
samd_ids <- read_csv("Gradu_AMR/samd_ids.csv")
samea_ids <- read_csv("Gradu_AMR/samea_ids.csv")
samn_ids <- read_csv("Gradu_AMR/samn_ids.csv")
srr_ids <- read_csv("Gradu_AMR/srr_ids.csv")
srs_ids <- read_csv("Gradu_AMR/srs_ids.csv")
srx_ids <- read_csv("Gradu_AMR/srx_ids.csv")
```
---

## 5. Match Metalog Sample IDs to TSE Object

```r
matches_err <- err_ids$ena_ers_sample_id[err_ids$ena_ers_sample_id %in% df_col$acc]
matches_ers <- ers_ids$ena_ers_sample_id[ers_ids$ena_ers_sample_id %in% df_col$acc]
matches_samd <- samd_ids$ena_ers_sample_id[samd_ids$ena_ers_sample_id %in% df_col$acc]
matches_samea <- samea_ids$ena_ers_sample_id[samea_ids$ena_ers_sample_id %in% df_col$acc]
matches_samn <- samn_ids$ena_ers_sample_id[samn_ids$ena_ers_sample_id %in% df_col$acc]
matches_srr <- srr_ids$ena_ers_sample_id[srr_ids$ena_ers_sample_id %in% df_col$acc]
matches_srs <- srs_ids$ena_ers_sample_id[srs_ids$ena_ers_sample_id %in% df_col$acc]
matches_srx <- srx_ids$ena_ers_sample_id[srx_ids$ena_ers_sample_id %in% df_col$acc]
```
Result: No matches found.

---

## 6. Try Matching by Numeric Suffix Only

```r
df_filtered <- df_col$acc[grepl("^(SRR|ERR|DRR)", df_col$acc)]

err_nums <- sub("^ERR", "", err_ids$ena_ers_sample_id)
srr_nums <- sub("^SRR", "", srr_ids$ena_ers_sample_id)
df_nums <- sub("^[A-Z]+", "", df_filtered)

err_df <- data.frame(type="ERR", id=err_ids$ena_ers_sample_id, num=err_nums)
srr_df <- data.frame(type="SRR", id=srr_ids$ena_ers_sample_id, num=srr_nums)
df_ref <- data.frame(type="df_col", id=df_filtered, num=df_nums)

combined <- rbind(err_df, srr_df)
matches <- merge(combined, df_ref, by="num", suffixes=c("_query","_df"))
matches
```
Result: No matches — confirmed manually.

---

## 7. Compare with SRA Metadata (Biosamples)

```r
SRA_metadata <- read.csv("~/F_AMR_project/Gradu_AMR/SRA_metadata_with_biosample.txt")
View(SRA_metadata)

# Check matches by biosample
matches_samd <- samd_ids$ena_ers_sample_id[samd_ids$ena_ers_sample_id %in% SRA_metadata$biosample]
matches_samea <- samea_ids$ena_ers_sample_id[samea_ids$ena_ers_sample_id %in% SRA_metadata$biosample]
matches_samn <- samn_ids$ena_ers_sample_id[samn_ids$ena_ers_sample_id %in% SRA_metadata$biosample]

# Match counts
length(matches_samd)  # 27
length(matches_samea) # 1617
length(matches_samn)  # 1976

# Total unique matches
all_matches <- c(matches_samd, matches_samea, matches_samn)
```
---

## 8. Annotate SRA Table with Metalog Matches

```r
SRA_metadata$Metalog <- ifelse(SRA_metadata$biosample %in% all_matches, 1, 0)

# Check duplicates
matched_rows <- SRA_metadata[SRA_metadata$Metalog == 1, ]
biosample_counts <- table(matched_rows$biosample)
duplicates <- names(biosample_counts[biosample_counts > 1])
length(duplicates)  # 381 duplicates

# Example duplicate
SRA_metadata %>% filter(biosample == "SAMEA2466887")
```
Observation: Same biosample can appear multiple times with different accession numbers — likely multiple samples from the same patient.

---

## 9. Filter Matched Samples

```r
SRA_metadata_matched <- SRA_metadata %>% filter(Metalog == 1)
View(SRA_metadata_matched)
```
---

## 10. Audit log

| Step | Input                             | Action                                | Output                     | Notes / Observations                                                   |
| ---- | --------------------------------- | ------------------------------------- | -------------------------- | ---------------------------------------------------------------------- |
| 1    | `TSE_filtered.rds`                | Load TSE object, extract column `acc` | `df_col`                   | Accession numbers extracted for matching                               |
| 2    | `human_sample_list.csv`           | Inspect large CSV with Unix           | `unique_sample_names.txt`  | Sample prefixes identified: ERR, ERS, SAMD, SAMEA, SAMN, SRR, SRS, SRX |
| 3    | `human_sample_list.csv`           | Split by prefix                       | `*_ids.csv`                | Easier to process smaller files in R                                   |
| 4    | `*_ids.csv`                       | Load into R                           | `err_ids`, `ers_ids`, etc. | Ready for matching                                                     |
| 5    | TSE acc vs Metalog IDs            | Direct match                          | `matches_*`                | **No matches found**                                                   |
| 6    | Numeric suffix match              | Match ignoring prefixes               | `matches`                  | **No matches**                                                         |
| 7    | `SRA_metadata_with_biosample.txt` | Match Metalog biosample IDs           | `all_matches`              | Found 3620 unique matches                                              |
| 8    | Annotate SRA metadata             | Add Metalog column                    | `SRA_metadata$Metalog`     | 5391 rows marked due to duplicates; 381 biosamples duplicated          |
| 9    | Filter matched rows               | Keep only Metalog = 1                 | `SRA_metadata_matched`     | Ready for downstream analysis                                          |
---

Conclusion:

* Original TSE object does not match Metalog IDs directly.
* Matching by numeric suffix also fails.
* Comparison with updated SRA biosample metadata yields 3620 unique matches, with some duplicates.
* Despite low coverage, keep this dataset for record and potential future use.

