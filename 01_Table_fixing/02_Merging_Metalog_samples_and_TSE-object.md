## 1.  Inspect SRA_metadata_with_biosample_corrected.txt


### 1.1 Examine Column Structure
```bash
head -n 1 SRA_metadata_with_biosample_corrected.txt | tr ',' '\n' | nl
```
Columns present:

1. acc
2. biosample
3. geo_loc_name_country_calc
4. geo_loc_name_country_continent_calc
5. platform
6. instrument
7. bioproject
8. avgspotlen
9. mbases
10. collection_date_sam

### 1.2 Extract Unique BioSample IDs from SRA Metadata
```bash
grep -oE '\bSAM(N|D|EA)[0-9]+' SRA_metadata_with_biosample_corrected.txt | sort -u > sra_biosample_ids.txt

wc -l sra_biosample_ids.txt
# count: 54178 
```

## 2. Inspect human_all_wide_2025-10-19.tsv.gz

This file contains wide-format SPIRE human metadata.

### 2.1 File Dimensions
```bash
gzcat human_all_wide_2025-10-19.tsv.gz | wc -l
# rows 85471
gzcat human_all_wide_2025-10-19.tsv.gz | head -n 1 | awk -F'\t' '{print NF}'
# columns: 2836
```
### 2.2 Confirm Presence of BioSample IDs

```bash
gzcat human_all_wide_2025-10-19.tsv.gz | grep -E '\bSAM(N|D|EA)[0-9]+' | wc -l
# 81436 rows contain BioSample IDs
```
### 2.3 Extract Unique BioSample IDs

```bash
gzcat human_all_wide_2025-10-19.tsv.gz | grep -oE '\bSAM(N|D|EA)[0-9]+' | sort -u > biosample_ids.txt
```
## 3. Identify Shared BioSample IDs

Find BioSample IDs present in both SRA and SPIRE datasets.
```bash

sort biosample_ids.txt > biosample_ids_sorted.txt
sort sra_biosample_ids.txt > sra_biosample_ids_sorted.txt

comm -12 biosample_ids_sorted.txt sra_biosample_ids_sorted.txt > matched_biosamples.txt

wc -l matched_biosamples.txt
# matches: 20339
```
## 4. Merge SRA and SPIRE Metadata

Overview of Strategy:
* Subset both datasets using shared BioSample IDs
* Merge on BioSample ID
* Retain clinically relevant metadata fields
* Remove irrelevant or sensitive columns
* Handle large files via chunked reading

### 4.1 Python Workflow
```python
import pandas as pd
import re
```
### 4.2 Define Input Files and Parameters

```python
human_file = "human_all_wide_2025-10-19.tsv.gz"
sra_file = "SRA_metadata_with_biosample_corrected.txt"
matched_file = "matched_biosamples.txt"
final_output = "cleaned_merged_final.tsv"

id_col_human = "spire_sample_name"
id_col_sra = "biosample"

chunksize = 20000  # safe for large files
```
### 4.3 Define Keywords for Column Filtering

These keywords retain columns related to:
* Demographics (age, sex, BMI)
* Diseases and conditions
* Antibiotic and drug exposure

```python
keywords = [
    "bmi", "body mass", "antibiotic", "antimicrobial", "drug",
    "sex", "gender", "age", "disease", "diagnosis", "condition",
    "infection", "immune", "inflammation", "therapy", "treatment",
    "cycline", "cillin", "phenicol", "bactam", "beta", "lactamase", "inhibitor",
    "cef", "loracarbef", "flomoxef", "latamoxef", "onam", "penem", "cilastatin",
    "cephalexin", "dazole", "mycin", "prim", "sulfa", "tylosin", "tilmicosin",
    "tylvacosin", "tildipirosin", "macrolide", "quinupristin", "dalfopristin",
    "xacin", "acid", "flumequine", "olaquindox", "sulfonamide", "tetracycline",
    "nitrofuran", "amphenicol", "polymyxin", "quinolone",
    "aminoglycoside", "teicoplanin", "vancin", "colistin", "polymyxin b",
    "nitro", "fura", "nifru", "micin", "tiamulin", "valnemulin",
    "xibornol", "clofoctol", "methenamine", "zolid", "lefamulins",
    "gepotidacin", "bacitracin", "novobiocin",
    "asthma", "cancer", "uti", "ibd", "crohn", "intestine",
    "ulcerative", "colitis"
]
keywords = [kw.lower() for kw in keywords]
```
### 4.4 Define Exclusion Keywords

Columns containing these terms are removed.

```python
exclude_keywords = [
    "sexual", "private", "identifier", "confidential", "contact",
    "gestational", "beverages", "stage", "vitamin", "alcohol",
    "average", "home", "sports", "strength", "siblings", "blockage"
]
```
### 4.5 Load Matched BioSample IDs

```python
matched_biosamples = pd.read_csv(
    matched_file, header=None, names=[id_col_sra], dtype=str
)
matched_set = set(matched_biosamples[id_col_sra])
print(f"Matched BioSamples: {len(matched_set)}")
# 20339
```
### 4.6 Subset SPIRE Human Metadata (Chunked)
```python

# 2. Subset human metadata by chunk
human_subset_file = "tmp_human_subset.tsv"
first = True

for chunk in pd.read_csv(human_file, sep="\t", dtype=str, chunksize=chunksize):
    filtered = chunk[chunk[id_col_human].isin(matched_set)]
    filtered.to_csv(
        human_subset_file,
        sep="\t",
        index=False,
        mode="w" if first else "a",
        header=first
    )
    first = False

human = pd.read_csv(human_subset_file, sep="\t", dtype=str)
print(f"Human subset: {human.shape}")
# (20339, 2836)
```
### 4.7 Subset SRA Metadata
```python
sra = pd.read_csv(sra_file, sep=",", dtype=str)
sra = sra[sra[id_col_sra].isin(matched_set)]

print(f"SRA subset: {sra.shape}")
# (24605, 10)
```
#### 4.8 Merge Datasets
```python
merged = pd.merge(
    human, sra,
    left_on=id_col_human,
    right_on=id_col_sra,
    how="inner",
    validate="many_to_many"
)
print(f"Merged: {merged.shape}")
# (24605, 2846)
```
Note:
More rows than unique BioSamples reflect multiple SRA accessions per BioSample.

### 4.9 Filter Columns by Relevance
Column selection rules
* Always keep SRA columns and BioSample ID
* Keep SPIRE columns containing relevant keywords
* Exclude columns with sensitive or irrelevant terms
```python
def keep_keyword(c):
    c_low = c.lower()
    return any(kw in c_low for kw in keywords)

def drop_excluded(c):
    c_low = c.lower()
    return any(bad in c_low for bad in exclude_keywords)

always_keep = set(sra.columns.tolist() + [id_col_human])
keyword_keep = {c for c in merged.columns if keep_keyword(c)}

cols_to_keep = always_keep | keyword_keep
cols_to_keep = sorted(
    c for c in cols_to_keep
    if c in merged.columns and not drop_excluded(c)
)

merged_clean = merged[cols_to_keep]
print(f"Final filtered columns: {len(cols_to_keep)}")
# 465
```
### 4.10 Save Final Output
```python
merged_clean.to_csv(final_output, sep="\t", index=False)
```

Final Output

cleaned_merged_final.tsv contains:
* All matched SRA accessions
* Linked SPIRE human metadata
* Demographic, clinical, and antibiotic-related variables
* Reduced from ~2800 to 465 curated columns

Summary:
* 20,339 shared BioSamples identified
* Many-to-many mapping preserved (multiple runs per sample)
* Metadata cleaned, filtered, and analysis-ready
