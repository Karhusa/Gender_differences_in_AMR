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

### 2.1 Python Script: unpack_jattr.py
```python
import pandas as pd
import json

# Read directly from the zipped file
df = pd.read_csv(
    "Sra_metadata_jun12_attributes.txt.zip",
    sep="\t",
    compression="zip",
    dtype={"jattr": "string"},
    low_memory=False
)

def parse_json_safe(s):
    try:
        return json.loads(s)
    except (json.JSONDecodeError, TypeError):
        return {}

jattr_dicts = df["jattr"].apply(parse_json_safe)

jattr_expanded = pd.json_normalize(jattr_dicts)

df_flat = pd.concat([df.drop(columns=["jattr"]), jattr_expanded], axis=1)

df_flat.to_csv("SRA_metadata_jun12_unpacked.csv", index=False)
```
### 2.2 Run the Script
```bash
python3 -m pip install pandas #version 2.3.3

python3 unpack_jattr.py
```

### 2.3 Inspect Output File
```bash

wc -l SRA_metadata_jun12.csv
# 17718 rows

awk -F',' '{print NF; exit}' SRA_metadata_jun12.csv
# 2132 columns

```
Interpretation:

* ~17k samples
* 2000 metadata fields due to expanded JSON attributes

## 3. Extract Gender and Sample Type

### 3.1 Load Metadata
```python

import pandas as pd
import ast

df = pd.read_csv("SRA_metadata_jun12.csv", low_memory=False)
```
## 3.2 Extract Gender Information
Identify candidate columns

```python
sex_cols = [c for c in df.columns if "sex" in c.lower() or "gender" in c.lower()]
print("Possible sex-related columns:", sex_cols)
```

Helper function to normalize gender values
```python
def extract_gender(val):
    if pd.isna(val):
        return None
    if isinstance(val, str) and val.strip().startswith("["):
        try:
            parsed = ast.literal_eval(val)
            if isinstance(parsed, list) and len(parsed) > 0:
                val = parsed[0]
        except:
            pass
    val = str(val).strip().lower()
    if "female" in val:
        return "female"
    if "male" in val:
        return "male"
    return None
```
Create a clean Gender column

```python
df["Gender"] = df[sex_cols].apply(lambda row: next(
    (extract_gender(v) for v in row if extract_gender(v) is not None), None), axis=1)

df = df.drop(columns=sex_cols)
```
### 3.3 Extract Sample Type (Fecal Samples)
Identify candidate columns

```python
sample_cols = [c for c in df.columns if "environment" in c.lower() or "sample" in c.lower()]
print("Possible sample-related columns:", sample_cols)
```
Helper function to detect fecal samples
```def detect_feces(val):
    if pd.isna(val):
        return None
    val = str(val).strip().lower()
    if any(keyword in val for keyword in ["stool", "feces", "human gut"]):
        return "feces"
    return None
```
Create Sample type column

```python
df["Sample type"] = df[sample_cols].apply(lambda row: next(
    (detect_feces(v) for v in row if detect_feces(v) is not None), None), axis=1)
```

### 3.4 Reorder Columns and Save Output

Move key annotation columns near the front.
```python
cols = list(df.columns)
for col in ["Gender", "Sample type"]:
    if col in cols:
        cols.insert(1, cols.pop(cols.index(col)))
df = df[cols]

df.to_csv("metadata_gender_sample.csv", index=False)
print("Created 'metadata_gender_sample.csv' with Gender and Sample type columns.")
```
metadata_gender_sample.csv contains:

* Original SRA metadata
* Expanded JSON attributes
* Clean Gender column (male, female, or NA)
* Clean Sample type column (feces or NA)

