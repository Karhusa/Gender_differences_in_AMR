Input files
* SRA_metadata_with_biosample_corrected (SRA biosample numbers from Katariina Pärnänen)
* SRA_metadata_gut_samples_full.csv (SRA Metadata)
* human_all_wide_2025-10-19.tsv.gz (Metalog metadata)

Objectives:
Combine SRA and Metalog metadata with common biosample numbers from SRA_metadata_with_biosample_corrected

___
## 1.  Merge SRA_metadata_with_biosample_corrected.txt

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


```python

import pandas as pd

corrected = pd.read_csv(
    "SRA_metadata_with_biosample_corrected.txt",
    sep=",",
    dtype=str
)

sra_full = pd.read_csv(
    "SRA_metadata_gut_samples_full.csv",
    sep=",",
    dtype=str
)

corrected.columns = corrected.columns.str.strip().str.lower()
sra_full.columns = sra_full.columns.str.strip().str.lower()

print("Corrected columns:", corrected.columns.tolist())
print("SRA full columns:", sra_full.columns.tolist())

if "acc" not in corrected.columns:
    raise ValueError("Column 'acc' not found in corrected file")

if "acc" not in sra_full.columns:
    raise ValueError("Column 'acc' not found in SRA full file")

print("Corrected acc duplicates:", corrected["acc"].duplicated().sum())
print("SRA full acc duplicates:", sra_full["acc"].duplicated().sum())

sra_combined = corrected.merge(
    sra_full,
    on="acc",
    how="left",
    suffixes=("", "_sra")
)

print("Rows in corrected:", len(corrected))
print("Rows after merge:", len(sra_combined))
print("Missing matches (acc not found in SRA full):", sra_combined["acc"].isna().sum())

sra_combined.to_csv(
    "SRA_combined_full.tsv",
    sep="\t",
    index=False
)

print("SRA_combined_full.tsv saved successfully.")
```
## Merge human_all_wide_2025-10-19.tsv.gz to SRA_combined_full.tsv

```python

metalog = pd.read_csv(
    "human_all_wide_2025-10-19.tsv.gz",
    sep="\t",
    compression="gzip",
    dtype=str
)

sra_combined = pd.read_csv(
    "SRA_combined_full.tsv",
    sep="\t",
    dtype=str
)

metalog.columns = metalog.columns.str.strip().str.lower()
sra_combined.columns = sra_combined.columns.str.strip().str.lower()

final_metadata = metalog.merge(
    sra_combined,
    left_on="spire_sample_name",
    right_on="biosample",
    how="inner",  # only keep rows with a match
    suffixes=("", "_sra")
)

final_metadata = final_metadata.dropna(axis=1, how='all')

cols = final_metadata.columns.tolist()
if "acc" in cols:
    cols.insert(0, cols.pop(cols.index("acc")))
final_metadata = final_metadata[cols]

final_metadata.to_csv(
    "Metalog_SRA_final_clean.tsv",
    sep="\t",
    index=False
)

num_rows, num_columns = final_metadata.shape
print(f"Final table saved: Metalog_SRA_final_clean.tsv")
print(f"Number of rows: {num_rows}")
print(f"Number of columns: {num_columns}")
```

Final table saved: Metalog_SRA_final_clean.tsv
Number of rows: 24605
Number of columns: 1378


### 4.2 Define Keywords for Column Filtering

These keywords retain columns related to:
* Demographics (age, sex, BMI)
* Diseases and conditions
* Antibiotic and drug exposure

```python
keywords = [
    "bmi", "body mass", "body", "antibiotic", "abx", "antimicrobial", "drug",
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
    "ulcerative", "colitis", "sx", "gndr", "female", "male", "uro", "dad", "mother",
    "mum", "maternal", "gestational", "sibling", "family"]


keywords = [kw.lower() for kw in keywords]
```
### 4.4 Define Exclusion Keywords

Columns containing these terms are removed.

```python
exclude_keywords = [
    "longitude", "latitude", "environment", "timepoint", "fmt", "hip", "waist", "note", "birth", "milk",
    "breakfast", "insulin", "measels", "physical", "gained", "diastolic", "drinks", "pizza",
    "comment", "insulin", "clotting", "blockage", "prognosis", "c-peptide", "constipation",
    "hemoglobin", "cholesterol", "arm"


]

study_code, sample_alias,collection_date, last_change, tax_id, weight_kg, height_cm, artificial, subject_id, intervention, available_info, pmid, doi, timeseries_available, timeseries_count, timeseries_duration, medication_with_parents, recipient_donor, birth_country, sample_description, birth_weight_kg, hip_cm, waist_cm, study_accession, bristol_stool_scale, "sample_title",  raw_metadata_age_at_diagnosis, animal_contact, raw_metadata_cr_(umol/l), raw_metadata_gi_malignancy_polyposis


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
### 4.6 Subset Metalog Human Metadata (Chunked)
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
* Keep Metalog columns containing relevant keywords
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
* Linked Metalog human metadata
* Demographic, clinical, and antibiotic-related variables
* Reduced from ~2800 to 465 curated columns

Summary:
* 20,339 shared BioSamples identified
* Many-to-many mapping preserved (multiple runs per sample)
* Metadata cleaned but remaining columns need to be inspected
