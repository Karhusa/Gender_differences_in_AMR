Input files:
* Sra_metadata_jun12_attributes.txt.zip

## 1. Inspect the Raw Metadata File

### 1.1 Preview File Contents and Filter Relevant Rows

The metadata file is compressed and very large, so we inspect it using Unix tools.

```bash
unzip -p Sra_metadata_jun12_attributes.txt.zip Sra_metadata_jun12_attributes.txt | head -n 5

```
### 1.2 extract all the column names:


```bash
unzip -p Sra_metadata_jun12_attributes.txt.zip | \
head -n 1 | \
tr '\t' '\n' | \
nl > column_names_numbered.txt
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
Interpretation:
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

gender_keywords = ["sex", "gender", "sx", "gndr", "female", "male", "mother", "maternal"]

gender_cols = [c for c in df.columns if any(k in c.lower() for k in gender_keywords)]
print("Possible sex/gender columns:", gender_cols)

df_gender_info = df[gender_cols]

df_gender_info.to_csv("SRA_sex_gender_columns.csv", index=False)

```
* Output files: columns_with_indexes.txt
* SRA_sex_gender_columns.csv
