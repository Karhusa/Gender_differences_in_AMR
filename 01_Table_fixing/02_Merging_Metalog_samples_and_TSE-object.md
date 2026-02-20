## 02_Merging_Metalog_samples_and_TSE-object.md

**Objectives:**
* Combine SRA and Metalog metadata with common biosample numbers from SRA_metadata_with_biosample_corrected
* Remove columns than include only "NA" values

Input files:
* SRA_metadata_with_biosample_corrected (SRA biosample numbers from Katariina Pärnänen)
* SRA_metadata_gut_samples_full.csv (SRA Metadata)
* human_all_wide_2025-10-19.tsv.gz (Metalog metadata)

___
## 1.  Merge SRA_metadata_with_biosample_corrected.txt and SRA_metadata_gut_samples_full.csv
* Both datafiles have column named "acc"

```python
import pandas as pd

corrected = pd.read_csv("SRA_metadata_with_biosample_corrected.txt",sep=",",dtype=st)

sra_full = pd.read_csv("SRA_metadata_gut_samples_full.csv",sep=",", dtype=str)

corrected.columns = corrected.columns.str.strip().str.lower()
sra_full.columns = sra_full.columns.str.strip().str.lower()

if "acc" not in corrected.columns:
    raise ValueError("Column 'acc' not found in corrected file")

if "acc" not in sra_full.columns:
    raise ValueError("Column 'acc' not found in SRA full file")

print("Corrected acc duplicates:", corrected["acc"].duplicated().sum())
print("SRA full acc duplicates:", sra_full["acc"].duplicated().sum())

sra_combined = corrected.merge(sra_full, on="acc", how="left", suffixes=("", "_sra"))

print("Rows in corrected:", len(corrected))
print("Rows after merge:", len(sra_combined))
print("Missing matches (acc not found in SRA full):", sra_combined["acc"].isna().sum())

sra_combined.to_csv("SRA_combined_full.tsv",sep="\t", index=False)

print("SRA_combined_full.tsv saved successfully.")
```
--

## 2. Merge human_all_wide_2025-10-19.tsv.gz to SRA_combined_full.tsv

### 2.1 Extract Metalog column names:
```bash
gzip -dc human_all_wide_2025-10-19.tsv.gz | \
head -n 1 | \
tr '\t' '\n' | \
nl -w2 -s$'\t' \
> Metalog_columns_with_indexes.txt
```
### 2.2 Merge datasets
```python

metalog = pd.read_csv("human_all_wide_2025-10-19.tsv.gz", sep="\t", compression="gzip", dtype=str)

sra_combined = pd.read_csv("SRA_combined_full.tsv",sep="\t",dtype=str)

metalog.columns = metalog.columns.str.strip().str.lower()
sra_combined.columns = sra_combined.columns.str.strip().str.lower()

final_metadata = metalog.merge(
    sra_combined,
    left_on="spire_sample_name",
    right_on="biosample",
    how="inner",  # only keep rows with a match
    suffixes=("", "_sra")
)

#remove columns thatn have only NA-values
final_metadata = final_metadata.dropna(axis=1, how='all')

# Move column named "acc" to index 0
cols = final_metadata.columns.tolist()
if "acc" in cols:
    cols.insert(0, cols.pop(cols.index("acc")))
final_metadata = final_metadata[cols]

final_metadata.to_csv("Metalog_SRA_final_clean.tsv",sep="\t", index=False)

num_rows, num_columns = final_metadata.shape
print(f"Final table saved: Metalog_SRA_final_clean.tsv")
print(f"Number of rows: {num_rows}")
print(f"Number of columns: {num_columns}")
```
Final table saved: Metalog_SRA_final_clean.tsv
Number of rows: 24605
Number of columns: 1378

## 3. Compress the final file

```bash
gzip Metalog_SRA_final_clean.tsv
```

Summary:
* 20,339 shared BioSamples identified
* Many-to-many mapping preserved (multiple runs per sample)
* Metadata cleaned but remaining columns need to be inspected
