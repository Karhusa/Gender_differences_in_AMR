# Metadata Cleaning Pipeline

## Overview

* Load and clean merged SRA metadata
* Remove empty / uninformative columns
* Variables were excluded if they were not relevant to the predefined research questions or planned statistical analyses, and therefore did not contribute meaningful information to the study objectives.
* Harmonise BMI, age, disease, infection, UTI, and antibiotic metadata
* Produce progressively cleaner TSV outputs (`kesken1.tsv`, `kesken2.tsv`)

Input files:
* Metalog_SRA_final_clean.tsv

---

## 1. Load Data

```python
import pandas as pd
import numpy as np
import re

df = pd.read_csv("Metalog_SRA_final_clean.tsv", sep="\t")
```

---

## 2. Remove empty columns

### 2.1 Columns With Na Values

```python
df = df.dropna(axis=1, how='all')
print(f" Shape: {df.shape}")
```
Shape: (24605, 1378)

### 2.2 Drop Columns Containing Only NaN or "No"

```python
def drop_nan_no_columns(df):
    cols_to_drop = []
    for col in df.columns:
        unique_vals = set(df[col].dropna().astype(str).str.strip())
        if len(unique_vals) == 0:
            cols_to_drop.append(col)
        elif unique_vals == {"No"}:
            cols_to_drop.append(col)
    df = df.drop(columns=cols_to_drop)
    return df, cols_to_drop

df, dropped = drop_nan_no_columns(df)
print(f" Shape: {df.shape}")
```
Shape: (24605, 1214)

---

## 3. Collect information about remaining columns

```python
summary_df = pd.DataFrame({
    "Column": df.columns,
    "Data_Type": df.dtypes,
    "Non_Null_Count": df.notnull().sum(),
    "Null_Count": df.isnull().sum(),
    "Unique_Values": df.nunique(),
})

summary_df.to_csv("Metadata_column_value_summary.tsv", sep="\t", index=False)
```
* No need to save, just for visual inspection

---

## 4. Filter columns

### 4.2 Define Keywords for Column Filtering

These keywords retain columns related to:
* Demographics (age, sex, BMI, family relationships)
* Geographic loaction
* Diseases and infections
* Antibiotic and drug exposure

```python
keywords = [
    "acc", "accession", "biosample", "bioproject",
    "age", "sex", "gender", "male", "female",
    "bmi", "body mass", "body",
    "region", "geographic", "location", "country", "contra", "urban", "rural"

    "family", "relationship", "sibling",
    "mother", "maternal", "mum", "mom", "dad", "father",
    "pregnant", "pregnancy", "gestational", "gravid",
    "menopau",

    "lifestyle", "life", "environmental",
    "diet", "smoke", "alcohol", "exercise",
    "disease", "diagnosis", "condition",
    "infection", "immune", "immunity", "inflammation",
    "cancer", "tumor", "malign", "malignancy",
    "asthma", "sepsis",
    "ibd", "crohn", "ulcerative", "colitis",
    "intestine", "gut",

    "uro", "urine", "tract", "uti",

    "therapy", "treatment", "medication", "drug",
    "antibiotic", "antibio", "abx", "antimicrobial",
    "vaccine", "vaccina",
    "biologics",

    "beta", "bactam", "lactamase", "inhibitor",
    "penem", "onam", "cilastatin",
    "cef", "cephalexin", "loracarbef", "flomoxef", "latamoxef",
    "cycline", "tetracycline",
    "cillin",
    "mycin", "micin", "macrolide",
    "aminoglycoside",
    "quinolone", "xacin", "flumequine",
    "sulfa", "sulfonamide",
    "phenicol", "amphenicol",
    "nitro", "nitrofuran", "fura", "nifru",
    "polymyxin", "polymyxin_b", "colistin",
    "vancin", "teicoplanin",
    "dazole",
    "prim",
    "tylosin", "tilmicosin", "tylvacosin", "tildipirosin",
    "quinupristin", "dalfopristin",
    "tiamulin", "valnemulin",
    "zolid", "lefamulins",
    "gepotidacin",
    "bacitracin", "novobiocin",
    "olaquindox",
    "xibornol", "clofoctol", "methenamine",

    "sx", "gndr"]

pattern = "|".join([re.escape(k) for k in keywords if k])

matched_columns = [
    col for col in df.columns
    if re.search(pattern, col, flags=re.IGNORECASE)]

extra_column = "raw_metadata_living_environment"
if extra_column in df.columns: matched_columns.append(extra_column)

matched_columns = list(set(matched_columns))
filtered_df = df[matched_columns]

print(f" Shape: {filtered_df.shape}")
```
Shape: (24605, 316)

* Output file: Filtered_Metadata.tsv (No need to save, just for inspection)


**Exclude columns related to:**
* technical metadata
* Lifestyle 
* Disease stage
* birth
* Unrelated treatments, diseases or measurements
* Too detailed geography/location

```python
exclude_keywords = [

    "16s", "metagenomic_sequencing", "nucleic_acid",
    "extractionprotocolid", "storageprotocol",
    "storage_temperature", "specimen_location_sam",
    "sample_location_full", "body_site_sam",
    "body_habitat_sam",

    "diet", "diet_sam", "diet_full", "raw_metadata_diet",
    "diet_as_a_baby", "particular_diet", "vegan",
    "lifestyle", "lifestylegeneral",
    "alcohol", "tobacco", "smoke", "smoked", "cigarette",
    "soft_drinks", "caffeinated", "physicalexercise", "excercise",

    "stage", "tnm_stage", "tumor_stage", "tumor_localization",
    "discovery_stage", "validation_stage", "hoehn_yahr_stage",
    "baseline_montreal_location", "disease_extent",
    "disease_activity", "disease_group", "disease_status_sam",

    "age_of_onset", "age_at_diagnosis",
    "diagnosis_date", "disease_duration",
    "disease_duration_years",

    "birth_gestational", "birth_country", "birth_outside",
    "day_of_life", "agedischargeddays",
    "feedingage", "neonate_pair",

    "treatment", "treatment_group", "treatment_duration",
    "treatment_effect", "treatment_where_sam",
    "treatment_by_sam", "treatment_(y_or_n)",
    "dosage", "dosage_1", "dosage_sam",
    "drug", "drug_statins", "drug_combine",
    "immunotherapy", "biologics",
    "inflammatory_drugs", "reaction_medication",
    "medication_code", "propranolol",

    "creatinine_mumol_l", "urea_nitrogen",
    "body_temp", "body_temperature",
    "body_%_fat", "body_fat_mass", "body_lean_mass",
    "age_z_score",

    "village", "village_sam", "village_grouping",
    "longitude", "country_birth",

    "submitter_id",

    "bilharzia", "schistosoma", "mumps_sam",

    "metadata_condition", "sexual_orientation",
    "infection_control_means"
]

pattern = "|".join([re.escape(k) for k in exclude_keywords])

cols_to_remove = filtered_df.columns[filtered_df.columns.str.contains(
    pattern, case=False, regex=True
)]

filtered_df = filtered_df.drop(columns=cols_to_remove)

print("Removed columns:", len(cols_to_remove))

print(f" Shape: {filtered_df.shape}")
```
Removed columns: 126
Shape: (24605, 190)

---
## 5.UTI-Related Columns

* Combine UTI-related columns into one UTI_History column by marking "Yes" if any column indicates a UTI, otherwise "No".
* Check results
* Remove unnecessary columns

```python
uti_cols = filtered_df.columns[filtered_df.columns.str.contains("uti", case=False)]
for col in uti_cols: print(f"{col}: {filtered_df[col].unique()}")

filtered_df["UTI_History"] = np.where(
    (filtered_df["raw_metadata_utis"].fillna(0) > 0) |
    (filtered_df["raw_metadata_history_of_recurrent_uti"].notna()) |
    (filtered_df["raw_metadata_ecoli_utis"].fillna(0) == 1) |
    (filtered_df["raw_metadata_diagnosed_utis"].fillna(0) == 1),
    "Yes", "No")

uri_cols = filtered_df.columns[filtered_df.columns.str.contains("urine", case=False)]
for col in uri_cols: print(f"{col}: {filtered_df[col].unique()}")

filtered_df.loc[filtered_df["raw_metadata_urineinfection"].fillna(0) == 1, "UTI_History"] = "Yes"

tract_cols = filtered_df.columns[filtered_df.columns.str.contains("tract", case=False)]
for col in tract_cols:print(f"{col}: {filtered_df[col].unique()}")

filtered_df.drop(columns=uti_cols, inplace=True)
filtered_df.drop(columns=uri_cols, inplace=True)
filtered_df.drop(columns=tract_cols, inplace=True)

filtered_df["UTI_History"].value_counts(dropna=False)
print(f" Shape: {filtered_df.shape}")
```
UTI_History
* No     24380
* Yes      225

Shape: (24605, 183)

---
## 6.Cancer-Related Columns

```python
cancer_cols = filtered_df.columns[filtered_df.columns.str.contains("cancer", case=False)]
for col in cancer_cols: print(f"{col}: {filtered_df[col].unique()}")

filtered_df.rename(columns={"raw_metadata_gi_cancer_past_3_months": "GI_Cancer"},inplace=True)
filtered_df["GI_Cancer"] = filtered_df["GI_Cancer"].notna().map({True: "Yes",False: "No"})

tumor_cols = filtered_df.columns[filtered_df.columns.str.contains("tumor", case=False)]
for col in tumor_cols: print(f"{col}: {filtered_df[col].unique()}")

filtered_df.rename(columns={"raw_metadata_tumor_location": "Colorectal_Cancer"}, inplace=True)
filtered_df["Colorectal_Cancer"] = filtered_df["Colorectal_Cancer"].notna().map({True: "Yes", False: "No"})

filtered_df["Colorectal_Cancer"].value_counts(dropna=False)
filtered_df["GI_Cancer"].value_counts(dropna=False)
```
Colorectal_Cancer
* No     24320
* Yes      285

GI_Cancer
* No     23938
* Yes      667
---

## 7. Infection-Related Columns

```python
infection_cols = filtered_df.columns[filtered_df.columns.str.contains("infection", case=False)]
for col in infection_cols: print(f"{col}: {filtered_df[col].unique()}")

filtered_df['Infection_history'] = pd.Series(dtype="object")
filtered_df['Infection'] = pd.Series(dtype="object")

infection_map = {
    'raw_metadata_trachealinfection': 'Tracheal Infection',
    'raw_metadata_bloodinfection': 'Blood Infection',
    'raw_metadata_gi_infection': 'GI Infection',
    'raw_metadata_otherinfection': 'Other Infection'}

def get_infections(row):
    infections = []
    for col, name in infection_map.items():
        if pd.notna(row[col]) and row[col] != 0:
            infections.append(name)
    return ', '.join(sorted(set(infections))) if infections else pd.NA

filtered_df['Infection'] = filtered_df.apply(get_infections, axis=1)

filtered_df.loc[filtered_df['Infection'].notna(), 'Infection_history'] = 'Yes'
filtered_df.loc[filtered_df['Infection'].isna(), 'Infection_history'] = 'No'

filtered_df["Infection_history"].value_counts(dropna=False)
filtered_df["Infection"].value_counts(dropna=False)

filtered_df.drop(columns=infection_cols, inplace=True)

print(f" Shape: {filtered_df.shape}")
```
Infection_history
* No     23663
* Yes      662
* NaN      280

Infection
* <NA>                                                    23943
* GI Infection                                              548
* Blood Infection                                            29
* Other Infection                                            20
* Blood Infection, Other Infection, Tracheal Infection       19
* Blood Infection, Tracheal Infection                        17
* Tracheal Infection                                         13
* Other Infection, Tracheal Infection                         9
* Blood Infection, Other Infection                            7

Shape: (24605, 181)

---

## 8.Location-Related Columns

```python
loc_cols = filtered_df.columns[filtered_df.columns.str.contains("location", case=False)]
for col in loc_cols: print(f"{col}: {filtered_df[col].unique()}")
```
* Many columns with mixed data which needs to be cleaned
* Clean and combine different location-related columns into under one column

```python
#Column name: geographic_location
filtered_df.rename(columns={"geographic_location": "Country"}, inplace=True)
filtered_df["Country"].value_counts(dropna=False)

# Column name: geographic_location__country_and_or_sea__sam
filtered_df["geographic_location__country_and_or_sea__sam"] = (
    filtered_df["geographic_location__country_and_or_sea__sam"]
    .str.strip("[]").str.replace("'", "", regex=False))

filtered_df["Country"] = filtered_df["Country"].fillna(filtered_df["geographic_location__country_and_or_sea__sam"])
filtered_df.drop(columns=["geographic_location__country_and_or_sea__sam"], inplace=True)

# Column name:raw_metadata_location
filtered_df["raw_metadata_location"] = filtered_df["raw_metadata_location"].replace("missing", pd.NA)
split_cols = filtered_df["raw_metadata_location"].str.split(":", n=1, expand=True)
filtered_df["Country"] = filtered_df["Country"].fillna(split_cols[0].str.strip())
filtered_df["City"] = split_cols[1].str.strip()
filtered_df["City"] = filtered_df["City"].replace({None: np.nan})

print(filtered_df["Country"].value_counts(dropna=False))
print("\n")
print(filtered_df["City"].value_counts(dropna=False))

filtered_df.drop(columns=["raw_metadata_location"], inplace=True)

# Column name: location
filtered_df["Continent"] = pd.NA 
filtered_df.loc[filtered_df["location"] == "North America", "Continent"] = "North America"
filtered_df.loc[filtered_df["location"] == "Judah", ["Country", "City"]] = pd.NA # what is "Judah"? lets just put it as NaN
filtered_df.loc[filtered_df["location"] == "Tanjung_Sepat", ["Country", "City"]] = ["Malaysia", "Tanjung_Sepat"]
mask = ~filtered_df["location"].isin(["North America", "Judah", "Tanjung_Sepat"]) & filtered_df["location"].notna()
split_loc = filtered_df.loc[mask, "location"].str.split(":", n=1, expand=True)

filtered_df.loc[mask, "Country"] = split_loc[0].str.strip()
filtered_df.loc[mask, "City"] = split_loc[1].str.strip()
filtered_df["City"] = filtered_df["City"].replace({None: np.nan})

print(filtered_df["Country"].value_counts(dropna=False))
print("\n")
print(filtered_df["City"].value_counts(dropna=False))

filtered_df.drop(columns=["location"], inplace=True)

# Column name: geographic_location__region_and_locality__sam
filtered_df["geographic_location__region_and_locality__sam"] = filtered_df[
    "geographic_location__region_and_locality__sam"].replace("nan", pd.NA)

filtered_df["geographic_location__region_and_locality__sam"] = (
    filtered_df["geographic_location__region_and_locality__sam"]
    .str.strip("[]").str.replace("'", "", regex=False))

split_loc = filtered_df["geographic_location__region_and_locality__sam"].str.split(":", n=1, expand=True)

filtered_df["Country"] = filtered_df["Country"].fillna(split_loc[0].str.strip())
filtered_df["City"] = filtered_df["City"].fillna(split_loc[1].str.strip())

print(filtered_df["Country"].value_counts(dropna=False))
print("\n")
print(filtered_df["City"].value_counts(dropna=False))

filtered_df.drop(columns=["geographic_location__region_and_locality__sam"], inplace=True)

# Column name: geographic_location__latitude__sam
filtered_df.drop(columns=["geographic_location__latitude__sam"], inplace=True)

cont_cols = filtered_df.columns[filtered_df.columns.str.contains("continent", case=False)]
for col in cont_cols: print(f"{col}: {filtered_df[col].unique()}")

# Clean and combine three continent related information columns
filtered_df["Continent"] = filtered_df["Continent"].astype("string")
for col in ["geo_loc_name_country_continent_calc_sra", "geo_loc_name_country_continent_calc"]:
    mask = filtered_df["Continent"].isna() | (filtered_df["Continent"] == "")
    filtered_df.loc[mask, "Continent"] = filtered_df.loc[mask, col]

filtered_df["Continent"] = filtered_df["Continent"].replace("uncalculated", pd.NA)
filtered_df["Continent"] = filtered_df["Continent"].str.title()
filtered_df.drop(columns=["geo_loc_name_country_continent_calc_sra", "geo_loc_name_country_continent_calc"], inplace=True)

filtered_df["Continent"].value_counts(dropna=False)

print(filtered_df["Continent"].value_counts(dropna=False))
print("\n")
print(filtered_df["Country"].value_counts(dropna=False))
print("\n")
print(filtered_df["City"].value_counts(dropna=False))
```

Clean the columns

```python

filtered_df["Country"] = filtered_df["Country"].replace([None, "NaN", "nan"], pd.NA)
filtered_df["City"] = filtered_df["City"].replace([None, "NaN", "nan"], pd.NA)

# Step 2: Define city -> country mapping
city_to_country = {
    "New York": "USA",
    "Bhopal, Madhya Pradesh": "India",
    "Nagpur": "India",
    "Kerala": "India",
    "Melghat": "India",
    "Kuala Lumpur": "Malaysia",
    "Hong Kong": "China",
    "Legong": "Indonesia",
    "Rasau": "Brunei"}

for city, country in city_to_country.items():
    mask = filtered_df["Country"] == city
    filtered_df.loc[mask, "City"] = city
    filtered_df.loc[mask, "Country"] = country

filtered_df["City"] = filtered_df["City"].replace([None, "NaN", "nan"], pd.NA)
filtered_df["Country"] = filtered_df["Country"].replace([None, "NaN", "nan"], pd.NA)

# Check results
filtered_df[["Country", "City", "Continent"]].value_counts(dropna=False)

# Add continents to some countries
country_continent_map = {
    "Fiji": "Oceania",
    "Israel": "Asia",
    "USA": "North America",
    "China": "Asia",
    "Sweden": "Europe",
    "Italy": "Europe",
    "Japan": "Asia"}

filtered_df.loc[
    filtered_df['Continent'].isna() & filtered_df['Country'].isin(country_continent_map.keys()),'Continent'
] = filtered_df.loc[
    filtered_df['Continent'].isna() & filtered_df['Country'].isin(country_continent_map.keys()), 'Country'
].map(country_continent_map)

print(f" Shape: {filtered_df.shape}")
```

---

## 9. Family-related columns

```python
family_cols = filtered_df.columns[filtered_df.columns.str.contains("family", case=False)]
for col in family_cols:
    print(f"{col}: {filtered_df[col].unique()}")

filtered_df["Parent"] = np.select(
    [filtered_df["host_family_relationship_sam"].str.contains("mother", case=False, na=False),
    filtered_df["host_family_relationship_sam"].str.contains("father", case=False, na=False)],
    ["Mother", "Father"],
    default="None")

filtered_df.drop(columns=["host_family_relationship_sam",], inplace=True)

maternal_cols = filtered_df.columns[filtered_df.columns.str.contains("maternal", case=False)]
for col in maternal_cols:
    print(f"{col}: {filtered_df[col].unique()}")
filtered_df.drop(columns=maternal_cols, inplace=True)

mo_cols = filtered_df.columns[filtered_df.columns.str.contains("mother", case=False)]
for col in mo_cols:
    print(f"{col}: {filtered_df[col].unique()}")

filtered_df.loc[
    (filtered_df["raw_metadata_mother"].notna()) | 
    (filtered_df["raw_metadata_mother_infant"].str.lower() == "mother"),
    "Parent"
] = "Mother"
filtered_df.drop(columns=mo_cols, inplace=True)

fa_cols = filtered_df.columns[filtered_df.columns.str.contains("father", case=False)]
for col in fa_cols:
    print(f"{col}: {filtered_df[col].unique()}")

filtered_df.loc[
    (filtered_df["raw_metadata_father"].notna()),
    "Parent"
] = "Father"
filtered_df.drop(columns=fa_cols, inplace=True)

filtered_df["Parent"].value_counts(dropna=False)

print(f" Shape: {filtered_df.shape}")
```
Parent
* None      24398
* Mother      163
* Father       44

Shape: (24605, 169)

---

## 10. Antibiotic Usage Processing 

```python
abc_cols = filtered_df.columns[filtered_df.columns.str.contains("antibio", case=False)]
for col in abc_cols:print(f"{col}: {filtered_df[col].unique()}")
abx_cols = filtered_df.columns[filtered_df.columns.str.contains("abx", case=False)]
for col in abx_cols:print(f"{col}: {filtered_df[col].unique()}")

def antibiotic_status(row):
    yes = False
    no = False

    if row.get('antibio_flag_6mo_sam') == 1 or row.get('antibio_flag_1mo_sam') == 1:
        yes = True
    if row.get('antibio_flag_6mo_sam') == 0 or row.get('antibio_flag_1mo_sam') == 0:
        no = True

    if row.get('raw_metadata_antibiotics_current') == 'Y':
        yes = True
    elif row.get('raw_metadata_antibiotics_current') == 'N':
        no = True

    if str(row.get('antibiotics_sam')).lower() == 'yes':
        yes = True
    elif str(row.get('antibiotics_sam')).lower() == 'no':
        no = True

    if row.get('raw_metadata_antibiotics_past_3_months') == 'Y':
        yes = True
    elif row.get('raw_metadata_antibiotics_past_3_months') == 'N':
        no = True

    if row.get('antibiotics_past_3months_sam') == 'Yes':
        yes = True
    elif row.get('antibiotics_past_3months_sam') in ['No', 'Not sure']:
        no = True

    val = row.get('antibiotics_last_6_months__y_n__sam')
    if isinstance(val, str):
        if val.lower() == 'yes':
            yes = True
        elif val.lower() == 'no':
            no = True

    val = row.get('raw_metadata_antibiotics_last3months')
    if isinstance(val, str):
        if 'yes' in val.lower():
            yes = True
        elif val.lower() == 'no':
            no = True

    if pd.notna(row.get('days_since_antibiotics')): yes = True
    if pd.notna(row.get('range_days_since_antibiotics')): yes = True

    if row.get('host_subject_assigned_to_antibiotics_arm_sam') == 'ABX+': yes = True
    elif row.get('host_subject_assigned_to_antibiotics_arm_sam') == 'ABX-': no = True

    if row.get('raw_metadata_antibiotics_at_birth') == 'YES': yes = True
    elif row.get('raw_metadata_antibiotics_at_birth') == 'NO': no = True

   
    if pd.notna(row.get('raw_metadata_antibiotics_with_admission_days')) and row.get('raw_metadata_antibiotics_with_admission_days') > 0:
        yes = True

    if pd.notna(row.get('raw_metadata_total_antibiotic_days')) and row.get('raw_metadata_total_antibiotic_days') > 0:
        yes = True
    elif row.get('raw_metadata_total_antibiotic_days') == 0:
        no = True

    numeric_abx_cols = [
        'raw_metadata_w_allabx',
        'raw_metadata_m_allabx',
        'raw_metadata_c_allabx',
        'raw_metadata_days_on_abx_14']

    for col in numeric_abx_cols:
        if pd.notna(row.get(col)) and row[col] > 0: yes = True
        elif row.get(col) == 0: no = True

    if yes:
        return "Yes"
    if no:
        return "No"

    return np.nan

filtered_df["Antibiotic_Use"] = filtered_df.apply(antibiotic_status, axis=1)

filtered_df.drop(columns=abx_cols, inplace=True)
filtered_df.drop(columns=abc_cols, inplace=True)

filtered_df["Antibiotic_Use"].value_counts(dropna=False)

print(f" Shape: {filtered_df.shape}")
```
Antibiotic_Use
* NaN    22388
* No      1273
* Yes      944

Shape: (24605, 151)



### 7.1 Identify Antibiotic Columns

```python
#Get all the antibiotic columns
antibiotic_cols_w = [c for c in filtered_df.columns if c.startswith("raw_metadata_w_")]
antibiotic_cols_c = [c for c in filtered_df.columns if c.startswith("raw_metadata_c_")]
antibiotic_cols_m = [c for c in filtered_df.columns if c.startswith("raw_metadata_m_")]

all_antibiotic_cols = antibiotic_cols_w + antibiotic_cols_c + antibiotic_cols_m

antibiotics = [
    "clindamycin", "amoxacillin", "cefazolin", "allabx", "vancomycin", "cefotaxime",
    "oxacillin", "ticarcillinclavulanate", "metronidazole", "penicilling",
    "ampicillinsulbactam", "gentamicin", "meropenem", "cefepime", "ampicillin",
    "cefoxitin", "sulfamethoxazoletrimethoprim"
]

for ab in antibiotics:
    w_col = f"raw_metadata_w_{ab}"
    c_col = f"raw_metadata_c_{ab}"
    m_col = f"raw_metadata_m_{ab}"
    new_col = ab.capitalize()

    def check_use(row):
        for col, cond in [(w_col, lambda x: x == 1),
                          (c_col, lambda x: x > 0),
                          (m_col, lambda x: x == 1)]:
            if col in filtered_df.columns:
                val = row[col]
                if pd.notna(val) and cond(val):
                    return "Yes"
        return "No"

    filtered_df[new_col] = filtered_df.apply(check_use, axis=1)


antibiotic_cols_cap = [ab.capitalize() for ab in antibiotics]

def update_antibiotics_used(row):
    if any(row[col] == 'Yes' for col in antibiotic_cols_cap if col in row):
        return 'Yes'
    return row['Antibiotic_Use']  # keep current value ('No')

filtered_df['Antibiotic_Use'] = filtered_df.apply(update_antibiotics_used, axis=1)

filtered_df = filtered_df.drop(columns=all_antibiotic_cols)


filtered_df["Antibiotic_Use"].value_counts(dropna=False)
print(f" Shape: {filtered_df.shape}")
```
Antibiotic_Use
NaN    22144
* No      1273
* Yes     1188

Shape: (24605, 120)

### 7.1 Macrolides

```python
filtered_df["received_macrolides_sam"].value_counts(dropna=False)
filtered_df.loc[filtered_df['received_macrolides_sam'] == 'Yes', 'Antibiotic_Use'] = "Yes"
filtered_df.rename(columns={'received_macrolides_sam': 'Macrolides'}, inplace=True)
filtered_df["Antibiotic_Use"].value_counts(dropna=False)

print(f" Shape: {filtered_df.shape}")
```
Antibiotic_Use
* NaN    22045
* Yes     1287
* No      1273


---
## 11 BMI Related columns

```python
bmi_cols = filtered_df.columns[filtered_df.columns.str.contains("bmi", case=False)]
for col in bmi_cols: print(f"{col}: {filtered_df[col].unique()}")
bo_cols = filtered_df.columns[filtered_df.columns.str.contains("body", case=False)]
for col in bo_cols: print(f"{col}: {filtered_df[col].unique()}")

# Combine all of the numeric bmi columns into one

def clean_bmi_column(col):
    col = col.astype(str).str.replace(',', '.', regex=False)
    col = pd.to_numeric(col, errors='coerce')
    col = col.where(col >= 10, np.nan)
    return col

filtered_df['bmi_clean'] = clean_bmi_column(filtered_df['bmi'])
filtered_df['bmi_sam_clean'] = clean_bmi_column(filtered_df['bmi_sam'])
filtered_df['host_bmi_clean'] = clean_bmi_column(filtered_df['host_body_mass_index_sam_s_dpl230'])
filtered_df['bmi_dpl92_clean'] = clean_bmi_column(filtered_df['body_mass_index_sam_s_dpl92'])

filtered_df['bmi_numeric'] = filtered_df[['bmi_clean', 'bmi_sam_clean', 
                                           'host_bmi_clean', 'bmi_dpl92_clean']].bfill(axis=1).iloc[:, 0]

def bmi_category(bmi):
    if pd.isna(bmi):
        return np.nan
    elif bmi < 18.5:
        return "Underweight (<18.5)"
    elif 18.5 <= bmi < 25:
        return "Normal (18.5-25)"
    elif 25 <= bmi < 30:
        return "Overweight (25-30)"
    else:
        return "Obese (>30)"

filtered_df['bmi_category'] = filtered_df['bmi_numeric'].apply(bmi_category)

filtered_df.loc[filtered_df['bmi_range'] == '>30', 'bmi_category'] = "Obese (>30)"
filtered_df.loc[filtered_df['bmi_range.2'] == '>30', 'bmi_category'] = "Obese (>30)"
filtered_df.loc[filtered_df['bmi_range.1'] == '20-25', 'bmi_category'] = "Normal (18.5-25)"
filtered_df.loc[filtered_df['bmi_range.1'] == '25-30', 'bmi_category'] = "Overweight (25-30)"
filtered_df.loc[filtered_df['bmi_range.3'] == '20-25', 'bmi_category'] = "Normal (18.5-25)"
filtered_df.loc[filtered_df['bmi_range.3'] == '25-30', 'bmi_category'] = "Overweight (25-30)"

cols_to_drop = [
    'bmi', 
    'bmi_sam', 
    'body_mass_index_sam_s_dpl92', 
    'host_body_mass_index_sam_s_dpl230',
    'bmi_range',
    'bmi_range.2',
    'bmi_range.1',
    'bmi_range.3'
]

filtered_df = filtered_df.drop(columns=cols_to_drop, errors='ignore')

filtered_df["bmi_category"].value_counts(dropna=False)
print(f" Shape: {filtered_df.shape}")
```
bmi_category
* NaN                    18324
* Normal (18.5-25)        3221
* Overweight (25-30)      1399
* Obese (>30)             1243
* Underweight (<18.5)      418

Shape: (24605, 118)


filtered_df = pd.read_csv("Filtered_Metadata.tsv", sep="\t")

## 12 IBD Related columns

```python
ibd_cols = filtered_df.columns[filtered_df.columns.str.contains("IBD", case=False)]
for col in ibd_cols: print(f"{col}: {filtered_df[col].unique()}")
cro_cols = filtered_df.columns[filtered_df.columns.str.contains("crohn", case=False)]
for col in cro_cols:print(f"{col}: {filtered_df[col].unique()}")
col_cols = filtered_df.columns[filtered_df.columns.str.contains("colitis", case=False)]
for col in col_cols: print(f"{col}: {filtered_df[col].unique()}")

filtered_df = filtered_df.rename(columns={"raw_metadata_ibd": "IBD"})

filtered_df.loc[
    (filtered_df['raw_metadata_crohns_disease'] == 'Y') & 
    (filtered_df['IBD'] != 'Yes'), 'IBD'] = 'Yes'

filtered_df['IBD'] = filtered_df['IBD'].replace({'Y': 'Yes', 'N': 'No'})

filtered_df['raw_metadata_colitis'] = filtered_df['raw_metadata_colitis'].replace({'Y': 'Yes', 'N': 'No'})
filtered_df['raw_metadata_necrotizingenterocolitis'] = filtered_df['raw_metadata_necrotizingenterocolitis'].replace({'Yes': 'Yes', 'No': 'No'})

for col in ['raw_metadata_colitis', 'raw_metadata_necrotizingenterocolitis']:
    filtered_df.loc[filtered_df['IBD'].isna() & filtered_df[col].notna(), 'IBD'] = \
        filtered_df[col]

filtered_df = filtered_df.drop(columns=['raw_metadata_colitis', 'raw_metadata_necrotizingenterocolitis', 'raw_metadata_crohns_disease'])

filtered_df['IBD'].value_counts(dropna=False)
print(f" Shape: {filtered_df.shape}")
```
IBD
* NaN    23523
* No      1024
* Yes       58

Shape: (24605, 115)

```python
columns_to_remove = [
    "env_package_sam",
    "raw_metadata_childhood_accomodation",
    "raw_metadata_accomodation",
    "raw_metadata_x3_in_the_past_2_weeks_have_you_used_an_oral_contrast",
    "best_response_to_therapy_sam",
    "raw_metadata_mumps",
    "raw_metadata_relationship_status",
    "raw_metadata_new_sexual_partner_last6months",
    "human_gut_environmental_package_sam",
    "raw_metadata_donor_relationship",
    "raw_metadata_type_of_exercise",
    "raw_metadata_average_time_per_time_(min)",
    "raw_metadata_body_fat_percentage",
    "raw_metadata_side_of_body_where_symptoms_first_appeared",
    "primary_search",
    "raw_metadata_bacillus_calmette_guerin_vaccine",
    "vaccine_name",
    "raw_metadata_oral_medication",
    "raw_metadata_vaccination_last3weeks",
    "raw_metadata_postmenstrualagedays",
    "raw_metadata_carbapenemase_bla_gene",
    "raw_metadata_totalantimicrobialsdays",
    "menopausal_status",
    "raw_metadata_housing_condition",
    "host_body_product_sam" 
]

filtered_df = filtered_df.drop(columns=columns_to_remove)
print(f" Shape: {filtered_df.shape}")
```
Shape: (24605, 90)


## 13. Pregnancy related columns

```python
preg_cols = filtered_df.columns[filtered_df.columns.str.contains("preg", case=False)]
for col in preg_cols: print(f"{col}: {filtered_df[col].unique()}")
gest_cols = filtered_df.columns[filtered_df.columns.str.contains("gest", case=False)]
for col in gest_cols: print(f"{col}: {filtered_df[col].unique()}")
grav_cols = filtered_df.columns[filtered_df.columns.str.contains("grav", case=False)]
for col in grav_cols: print(f"{col}: {filtered_df[col].unique()}")


filtered_df['Pregnant'] = np.nan

columns_to_binary = [
    'raw_metadata_gestational_age_weeks',
    'pregnancy_week',
    'gestational_age_weeks',
    'gestational_state',
    'raw_metadata_gravida']

filtered_df['Pregnant'] = np.where(
    filtered_df[columns_to_binary].notna().any(axis=1),
    'Yes', 'No')

filtered_df.drop(columns=columns_to_binary, inplace=True)


filtered_df['Pregnant'].value_counts(dropna=False)
print(f" Shape: {filtered_df.shape}")
```
Pregnant
* No     23125
* Yes     1480

Shape: (24605, 86)

```

13. Gender

```python
sex_cols = filtered_df.columns[filtered_df.columns.str.contains("sex", case=False)]
for col in sex_cols:print(f"{col}: {filtered_df[col].unique()}")
gen_cols = filtered_df.columns[filtered_df.columns.str.contains("gen", case=False)]
for col in gen_cols: print(f"{col}: {filtered_df[col].unique()}")

for col in ['host_sex_sam', 'gender_sam']:
    filtered_df[col] = filtered_df[col].str.replace(r"[\[\]']", "", regex=True)

filtered_df['Sex'] = np.where(
    filtered_df['sex'].notna(),
    filtered_df['sex'],
    np.where(
        filtered_df['sex_calc'].notna(),
        filtered_df['sex_calc'],
        np.where(
            filtered_df['host_sex_sam'].notna(),
            filtered_df['host_sex_sam'],
            np.nan)))

filtered_df['Sex'] = filtered_df['Sex'].replace({'not collected': np.nan, 'not provided': np.nan})

columns_to_drop = ['host_sex_sam', 'sex_calc', 'sex', 'gender_sam']
filtered_df.drop(columns=columns_to_drop, inplace=True)

filtered_df['Sex'].value_counts(dropna=False)
print(f" Shape: {filtered_df.shape}")

```
Sex
* NaN       9824
* female    7430
* male      7351

print(f" Shape: {filtered_df.shape}")
Shape: (24605, 83)


filtered_df.to_csv("Filtered_Metadata.tsv", sep="\t", index=False)
filtered_df = pd.read_csv("Filtered_Metadata.tsv", sep="\t")
