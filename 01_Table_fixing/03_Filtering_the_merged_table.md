# Metadata Cleaning Pipeline

## Overview

* Load and clean merged SRA metadata
* Remove empty / uninformative columns
* Variables were excluded if they were not relevant to the predefined research questions or planned statistical analyses, and therefore did not contribute meaningful information to the study objectives.
* Harmonise BMI, age, disease, infection, UTI, and antibiotic metadata
* Produce progressively cleaner TSV outputs (`kesken1.tsv`, `kesken2.tsv`)

Input files:
* Filtered_Metadata.tsv

---

## 1. Load Data

```python
import pandas as pd
import numpy as np
import re

df = pd.read_csv("Filtered_Metadata.tsv", sep="\t")
```

---

## 2. Remove empty columns.

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
print("Dropped columns:", dropped)
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
    "region", "geographic", "location", "country", "contra"

    "family", "relationship", "sibling",
    "mother", "maternal", "mum", "mom", "dad", "father",
    "pregnant", "pregnancy", "gestational", "gravid",
    "menopau",

    "lifestyle", "life", "environmental",
    "diet", "smoke", "alcohol", "exercise",
    "disease", "diagnosis", "condition",
    "infection", "immune", "immunity", "inflammation",
    "cancer", "tumor", "malign", "malignancy",
    "asthma",
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

    "sx", "gndr"
]

pattern = "|".join([re.escape(k) for k in keywords if k])

matched_columns = [
    col for col in df.columns
    if re.search(pattern, col, flags=re.IGNORECASE)
]

extra_column = "raw_metadata_living_environment"

if extra_column in df.columns:
    matched_columns.append(extra_column)

matched_columns = list(set(matched_columns))

filtered_df = df[matched_columns]

# filtered_df.to_csv("Filtered_Metadata.tsv", sep="\t", index=False)

print(f" Shape: {filtered_df.shape}")
```
Shape: (24605, 317)

* Output file: Filtered_Metadata.tsv
     * No need to save, just for inspection



Exclude columns related to:
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

cols_to_remove = df.columns[df.columns.str.contains(
    pattern, case=False, regex=True
)]

df = df.drop(columns=cols_to_remove)

print("Removed columns:", len(cols_to_remove))

print(f" Shape: {df.shape}")
 Shape: (24605, 191)
```

## 4. BMI Cleaning and Categorization

### 4.1 Inspect BMI Columns

```python
bmi_cols = df.columns[df.columns.str.contains("bmi", case=False)]

for col in bmi_cols:
    print(f"{col}: {df[col].unique()}")


```

### 4.2 Convert BMI to Numeric and Create Categories

* Numeric BMI column

```python

df["bmi"] = pd.to_numeric(df["bmi"], errors="coerce")
df["bmi_sam"] = pd.to_numeric(df["bmi_sam"], errors="coerce")

df["bmi"] = df["bmi"].round(1)
df["bmi_sam"] = df["bmi_sam"].round(1)

df["BMI"] = df["bmi"].fillna(df["bmi_sam"])

conflicts = df[
     df["bmi"].notna() &
     df["bmi_sam"].notna() &
     (df["bmi"] != df["bmi_sam"])
]
 
conflicts.shape

df["BMI"].count()
```
No conflicts

BMI has 5774 values



```python
df["bmi"] = pd.to_numeric(df["bmi"], errors="coerce")

bins = [0, 18.5, 25, 30, float("inf")]
labels = [
    "Underweight (<18.5)",
    "Normal (18.5-25)",
    "Overweight (25-30)",
    "Obese (>30)"
]

df["BMI_range_new"] = pd.cut(df["bmi"], bins=bins, labels=labels)


```

---

### 4.3 Override Using Text-Based `bmi_range`

```python
bmi_range_map = {
    "16-20": "Underweight (<18.5)",
    "18.5-28": "Normal (18.5-25)",
    "20-25": "Normal (18.5-25)",
    "25-30": "Overweight (25-30)",
    "over 25": "Overweight (25-30)",
}

mask = df["bmi_range"].notna()
df.loc[mask, "BMI_range_new"] = df.loc[mask, "bmi_range"].map(bmi_range_map)

filtered_df = df.drop(columns=cols_to_remove)
```

---

### 4.4 Add Additional Categories

```python
if not pd.api.types.is_categorical_dtype(df["BMI_range_new"]):
    df["BMI_range_new"] = df["BMI_range_new"].astype("category")

new_categories = ["Obese (>30)", "Normal/Overweight (<30)"]
existing = df["BMI_range_new"].cat.categories

to_add = [cat for cat in new_categories if cat not in existing]

df["BMI_range_new"] = df["BMI_range_new"].cat.add_categories(to_add)
```

---

### 4.5 Override Using `BMI_range`

```python
BMI_range_map = {">30": "Obese (>30)", "<30": "Normal/Overweight (<30)"}

mask2 = df["BMI_range"].notna()
df.loc[mask2, "BMI_range_new"] = df.loc[mask2, "BMI_range"].map(BMI_range_map)
```

---

### 4.6 Review and Cleanup

```python
print(df["BMI_range_new"].value_counts(dropna=False))
print(df[["bmi", "bmi_range", "BMI_range", "BMI_range_new"]].head(10))

df = df.drop(columns=bmi_cols)
```

---

## 5. Drop Additional Unnecessary Metadata Columns

```python
columns_to_drop = [
    "raw_metadata_BloodInfection",
    "raw_metadata_age-block",
    "raw_metadata_age_of_onset",
    "raw_metadata_body_fat_percentage",
    "raw_metadata_height_for_age_z_score",
    "raw_metadata_maternal_age_at_delivery_years",
    "raw_metadata_treatment_batch",
    "raw_metadata_PostmenstrualAgeDays",
    "raw_metadata_MaternalAgeYears",
    "raw_metadata_Metagenomic_sequencing________(Y_or_N)",
    "raw_metadata_age_group_16S",
    "raw_metadata_storageprotocol",
    "village",
    "raw_metadata_AgeDischargedDays",
    "raw_metadata_Age_at_collection",
    "raw_metadata_Age_at_diagnosis",
    "raw_metadata_Condition",
    "raw_metadata_Drug_insulin",
    "raw_metadata_Drug_statins",
    "raw_metadata_Treatment_duration_(months)",
    "raw_metadata_age_at_diagnosis",
    "raw_metadata_diagnosis_date",
    "raw_metadata_disease_duration",
    "raw_metadata_disease_duration_year",
    "raw_metadata_disease_duration_years",
    "raw_metadata_disease_extent",
    "raw_metadata_housing_condition",
    "raw_metadata_treatment_group",
    "raw_metadata_treatment_effect",
    "raw_metadata_treatment_batch",
    "raw_metadata_response_to_treatment",
    "raw_metadata_mental_illness_diagnosis",
    "raw_metadata_Blood_urea_nitrogen",
    "raw_metadata_Drug_propranolol",
    "raw_metadata_NecrotizingEnterocolitis",
    "raw_metadata_Treatment_(Y_or_N)",
    "raw_metadata_previous_bilharzia_treatment",
    "raw_metadata_Anti_inflammatory_drugs",
]

df = df.drop(columns=columns_to_drop)
print(f" Shape: {df.shape}")
```

---

## 6. Remove Acid-Related Columns

```python
acid_cols = df.columns[df.columns.str.contains("acid", case=False)]
# 'raw_metadata_Alpha_lipoc_acid', 'raw_metadata_Total_bile_acid', 
# 'raw_metadata_Uric_acid', 'raw_metadata_valproic_acid'], dtype='object'
df = df.drop(columns=acid_cols)

```

---

## 7.Cancer-Related Columns

```python
cancer_cols = df.columns[df.columns.str.contains("cancer", case=False)]
#raw_metadata_GI_cancer_past_3_months: [nan 'N' 'Y']
```

---

## 8. Infection-Related Columns
### 8.1 Inspecto Infection Columns 

```python
infection_cols = [c for c in df.columns if c.startswith("raw_metadata_Infection_")]

for col in infection_cols:
    print(f"{col}: {df[col].unique()}")
#raw_metadata_Infection_Control_Means: [nan 'Cohorted' 'Isolated' 'Outpatient']

df = df.drop(columns=infection_cols)

```

---

## 8.2 Inspect Other Infection Columns
```python

infection_cols1 = df.columns[df.columns.str.contains("infection", case=False)]
for col in infection_cols1:
    print(f"{col}: {df[col].unique()}")

#raw_metadata_GI_infection: [nan 'No' '3-4 years ago, stomach bug']
#raw_metadata_OtherInfection: [nan  0.  1.]
#raw_metadata_TrachealInfection: [nan  0.  1.]
#raw_metadata_UrineInfection: [nan  0.  1.]
#raw_metadata_multi_drug_resistant_organism_infection: [nan 'Negative' 'Positive']
#raw_metadata_schistosoma_infection_intensity: ...

### 8.3 Drop Selected Columns

df = df.drop(columns=[
    "raw_metadata_OtherInfection",
    "raw_metadata_TrachealInfection",
    "raw_metadata_schistosoma_infection_intensity",
])
```

---

# 9. Antibiotic Usage Processing

## 9. Antibiotic Usage Processing ### 9.1 Identify Antibiotic Columns

```python
antibiotic_cols_w = [c for c in df.columns if c.startswith("raw_metadata_w_")]
antibiotic_cols_c = [c for c in df.columns if c.startswith("raw_metadata_c_")]
antibiotic_cols_m = [c for c in df.columns if c.startswith("raw_metadata_m_")]

all_antibiotic_cols = antibiotic_cols_w + antibiotic_cols_c + antibiotic_cols_m
```

## 9.2 Create Antibiotic List per Sample

```python
def get_antibiotics_list(row):
    used = []
    for col in all_antibiotic_cols:
        if row[col] == 1:
            name = col.replace("raw_metadata_w_", "") \
                      .replace("raw_metadata_c_", "") \
                      .replace("raw_metadata_m_", "")
            used.append(name)
    return used 

df["antibiotics_list"] = df.apply(get_antibiotics_list, axis=1)
```
---

## 9.3 Expand Lists Into Columns

```python

antibiotics_expanded = pd.DataFrame(
    df["antibiotics_list"].to_list(),
    index=df.index
).rename(lambda x: f"antibiotic_{x+1}", axis=1)

df = df.join(antibiotics_expanded)
df["name_of_antibiotic"] = df.apply(get_antibiotics_list, axis=1)

# Sanity check
df["name_of_antibiotic"].apply(tuple).unique()
```

## 9.4 Add Yes/No Column wheter antibiotics are used or not

```python
df["Antibiotics_used"] = df["antibiotics_list"].apply(lambda x: "Yes" if len(x) > 0 else "No")
df["Antibiotics_used"].value_counts(dropna=False)
# Antibiotics_used
# No     24459
# Yes      146
```
---

## 9.5 Remove Raw Antibiotic Columns

```python
df = df.drop(columns=all_antibiotic_cols)
```
# 10. Disease Columns 

## 10.1 Inspect Disease Columns

```python
disease_cols = df.columns[df.columns.str.contains("disease", case=False)]
for col in disease_cols:
    print(f"{col}: {df[col].unique()}")

#raw_metadata_Celiac_disease: [nan 'No' 'N' 'Y']
#raw_metadata_Crohns_disease: [nan 'N' 'Y']
#raw_metadata_Disease_activity_(Y_or_N): [nan 'Y' 'N']
#raw_metadata_Intestinal_disease: [nan 'N' 'Y']
#raw_metadata_diagnosed_with_disease: [nan 'No' 'Yes' 'not provided']
#raw_metadata_disease_cause: [nan 'HBV,alcohol' 'HBV' 'alcohol' 'Hepatitis C virus related''HBV,Hepatitis E virus related']
#raw_metadata_disease_group: [nan 'Healthy' 'Stage_III_IV' 'Stage_I_II' 'MP' 'Stage_0' 'HS']
#raw_metadata_diseases: [nan 'NEC' 'cellulitis,MRSA' 'acinetobacter anitratus' 'sepsis''cellulitis' 'bradycardia']
#raw_metadata_host_disease: [nan 'Acute Lymphoblastic Leukemia' 'Acute Myeloid Leukemia']
#raw_metadata_subject_disease_status_full: [nan 'COHORT' 'gestational diabetes mellitus']
#subject_disease_status: ['CTR' 'atopy' 'COHORT' "Crohn's disease" 'type 2 diabetes' ... many more ]
#subject_disease_status_full: [nan 'preeclampsia' 'mild preeclampsia' 'oligohydramnios' ... many more]
```

## 10.2 Drop Selected Disease Columns

```python
columns_to_drop = [
    "raw_metadata_Celiac_disease",
    "raw_metadata_Disease_activity_(Y_or_N)",
    "raw_metadata_diagnosed_with_disease",
    "raw_metadata_disease_cause",
    "raw_metadata_disease_group"
]

df = df.drop(columns=columns_to_drop)
```
---

# 11. Remove Additional Columns

```python

df = df.drop(columns=["raw_metadata_weight_for_age_z_score"])
print(df.shape)
# (24605, 63)
```
---

# 12. IBD Status

```python
df["IBD"] = (
    df["raw_metadata_IBD"]
    .map({"Y": "Yes", "N": "No"})
    .astype("category")
)

df = df.drop(columns=["raw_metadata_IBD"])
```
---

# 13. Export Intermediate Table
```python
df.to_csv("kesken1.tsv", sep="\t", index=False)
```
---

# 14. UTI History

## 14.1 Inspect UTI related columns

```python
uti_cols = [col for col in df.columns if "uti" in col.lower()]

for col in uti_cols:
    print(f"Column: {col}")
    print(df[col].unique())
    print()
```
---

## 14.2 Create a single UTI history variable.
* We also have one column (raw_metadata_UrineInfection) from 8.2 Inspect Remaining Infection Columns
  
```python
df["UTI_history"] = np.where(
    (df.get("raw_metadata_UTIs", 0) > 0) |
    (df.get("raw_metadata_diagnosed_utis", 0) == 1) |
    (df.get("raw_metadata_ecoli_utis", 0) == 1) |
    (df.get("raw_metadata_history_of_recurrent_uti") == "Recurrent UTIs") |
    (df.get("raw_metadata_UrineInfection") == 1),
    "Yes", "No"
)

df = df.drop(columns=[
    "raw_metadata_UTIs",
    "raw_metadata_diagnosed_utis",
    "raw_metadata_ecoli_utis",
    "raw_metadata_history_of_recurrent_uti",
    "raw_metadata_UrineInfection",
])

print(df["UTI_history"].value_counts(dropna=False))


```
UTI_history
* No     24380
* Yes      225

---

#15. Final Antibiotic Harmonization

* Infer antibiotic exposure from multiple metadata sources.

```python

antibiotic_cols = df.columns[df.columns.str.contains("antibiotic", case=False)]

for col in antibiotic_cols:
    print(f"\n=== {col} ===")
    print(df[col].value_counts(dropna=False).head(10))

# look through the columns and remove unnecessary columns with non-important information
exclude_cols = ['antibiotics_list', 'name_of_antibiotic', 'Antibiotics_used'] + [f'antibiotic_{i}' for i in range(1, 10)]

filtered_antibiotic_cols = [col for col in antibiotic_cols if col not in exclude_cols]

for col in filtered_antibiotic_cols:
    print(f"\n=== {col} ===")
    print(df[col].value_counts(dropna=False).head(10))

# Drug_antibiotic_last3y
# days_since_antibiotics
# range_days_since_antibiotics
# raw_metadata_Antibiotics_Last3months
# raw_metadata_Antibiotics_current
# raw_metadata_Antibiotics_past_3_months
# raw_metadata_antibiotic_use
# raw_metadata_antibiotics_at_birth
# raw_metadata_antibiotics_with_admission_days
# raw_metadata_total_antibiotic_days


df.loc[df["Drug_antibiotic_last3y"] == "2 months", "Antibiotics_used"] = "Yes"
df.loc[df["days_since_antibiotics"].notna(), "Antibiotics_used"] = "Yes"
df.loc[df["range_days_since_antibiotics"].notna(), "Antibiotics_used"] = "Yes"
df.loc[df["raw_metadata_Antibiotics_Last3months"].str.contains("Yes", na=False), "Antibiotics_used"] = "Yes"
df.loc[df["raw_metadata_Antibiotics_current"] == "Y", "Antibiotics_used"] = "Yes"
df.loc[df["raw_metadata_Antibiotics_past_3_months"] == "Y", "Antibiotics_used"] = "Yes"
df.loc[df["raw_metadata_antibiotics_at_birth"].str.upper() == "YES", "Antibiotics_used"] = "Yes"
df.loc[df["raw_metadata_antibiotics_with_admission_days"].gt(0), "Antibiotics_used"] = "Yes"
df.loc[df["raw_metadata_total_antibiotic_days"].gt(0), "Antibiotics_used"] = "Yes"


# Remove intermediate antibiotic columns:

df = df.drop(columns=[
    "Drug_antibiotic_last3y",
    "days_since_antibiotics",
    "range_days_since_antibiotics",
    "raw_metadata_Antibiotics_Last3months",
    "raw_metadata_Antibiotics_current",
    "raw_metadata_Antibiotics_past_3_months",
    "raw_metadata_antibiotic_use",
    "raw_metadata_antibiotics_at_birth",
    "raw_metadata_antibiotics_with_admission_days",
    "raw_metadata_total_antibiotic_days",
])
```
---

# 17. Final Export

df.to_csv("kesken2.tsv", sep="\t", index=False)

---

## Outputs

* Metadata_column_value_summary.tsv
* `kesken2.tsv`

---

df["raw_metadata_sexual_orientation"].unique()

# Antibiotics: 1.1 S01AA








