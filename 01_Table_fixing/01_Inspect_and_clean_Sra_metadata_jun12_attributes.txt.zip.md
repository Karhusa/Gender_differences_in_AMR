## 01_Inspect_and_clean_Sra_metadata_jun12_attributes.txt.zip.md

**Objectives:**
* Filter age and gut related samples from the whole dataset.

Input files:
* Sra_metadata_jun12_attributes.txt.zip

## 1. Inspect the Raw Metadata File

### 1.1 Preview File Contents and Filter Relevant Rows

* The metadata file is compressed and very large, so we inspect it using Unix tools.

```bash
unzip -p Sra_metadata_jun12_attributes.txt.zip Sra_metadata_jun12_attributes.txt | head -n 5
```
### 1.2 Extract all the column names:

```bash
unzip -p Sra_metadata_jun12_attributes.txt.zip | \
head -n 1 | \
tr '\t' '\n' | \
nl -w2 -s$'\t' > SRA_columns_with_indexes.txt
```

## 2. Unpack the jattr Column (JSON Attributes)

The jattr column contains nested JSON-like metadata. We unpack it into individual columns using Python.

### 2.1 Unzip

```bash
unzip Sra_metadata_jun12_attributes.txt.zip
```

### 2.2 Python Script: unpack_jattr.py
```python
import pandas as pd
import json

txt_file = "Sra_metadata_jun12_attributes.txt"

df = pd.read_csv(txt_file, sep=",", dtype="string", low_memory=False)

df.columns = df.columns.str.strip()

if "jattr" not in df.columns:
    raise ValueError(f"Column 'jattr' not found! Available columns: {df.columns.tolist()}")

def parse_json_safe(s):
    try:
        return json.loads(s)
    except (json.JSONDecodeError, TypeError):
        return {}

jattr_dicts = df["jattr"].apply(parse_json_safe)

jattr_expanded = pd.json_normalize(jattr_dicts)

df_flat = pd.concat([df.drop(columns=["jattr"]), jattr_expanded], axis=1)

df_flat.to_csv("SRA_metadata_jun12_unpacked.csv", index=False)

print("Flattened CSV saved successfully. Rows:", len(df_flat))

```
### 2.3 Run the Script
```bash
source venv2/bin/activate

python3 -m pip install pandas

python3 unpack_jattr.py
```
Result:
* Output file: SRA_metadata_jun12_unpacked.csv

### 2.4 Inspect Output File
```bash
wc -l SRA_metadata_jun12_unpacked.csv
awk -F',' '{print NF; exit}' SRA_metadata_jun12_unpacked.csv
```
Results:
* 69916 rows
* 2395 columns

---

## 3. Extract Gender and Sample Type

## 3.1 Extract Gender Information

* Identify candidate columns
```python

df = pd.read_csv(
    "SRA_metadata_jun12_unpacked.csv",
    sep=",",
    dtype="string",
    low_memory=False
)

with open("columns_with_indexes.txt", "w") as f:
    for i, col in enumerate(df.columns):
        f.write(f"{i}: {col}\n")

sex_keywords = ["sex", "gender", "sx", "gndr", "female", "male", "mother", "maternal"]

sex_cols = [c for c in df.columns if any(k in c.lower() for k in sex_keywords)]
print("Sex/gender-related columns found:", sex_cols)

df_known_sex = df[df[sex_cols].notna().any(axis=1)]

df_known_sex_all = df_known_sex.copy()

df_known_sex_all.to_csv("SRA_metadata_known_sex_full.csv", index=False)

print(f"Saved {len(df_known_sex_all)} samples with known sex to 'SRA_metadata_known_sex_full.csv'")
```
* No need to save this file, because the gut samples need to be separates.
* This was for visual inspection.


## 3.2 Extract Sample type Information

```r
df = pd.read_csv( "SRA_metadata_known_sex_full.csv", sep=",", dtype="string", low_memory=False)

gut_keywords = [
    "stool", "feces", "faeces",
    "human gut", "human-gut", "human gut microbiome",
    "stool sample", "human stool",
    "gastric biopsy", "lower digestive tract",
    "appendix", "stomach", "antrum", "corpus",
    "ileum", "colon", "intestine"
]

gut_keywords = [k.lower() for k in gut_keywords]

mask = df.apply(
    lambda row: any(
        any(keyword in str(cell).lower() for keyword in gut_keywords)
        for cell in row),
    axis=1)

df_gut_samples = df[mask].copy()

df_gut_samples.to_csv("SRA_metadata_gut_sex_samples_full.csv", index=False)

print(f"Saved {len(df_gut_samples)} gut-related samples.")
```

**Output files: **
* SRA_columns_with_indexes.txt
* SRA_metadata_known_sex_full.csv
* df_gut_samples.to_csv("SRA_metadata_gut_sex_samples_full.csv", index=False)
    * file is too large for github 

