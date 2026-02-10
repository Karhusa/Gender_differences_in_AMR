# Headline

## 1. Size and Shape

```python
import numpy as np
import pandas as pd
import re

df = pd.read_csv("kesken4.tsv", sep="\t")

print("DataFrame shape:", df.shape)
# DataFrame shape: (24605, 39)

print(df.info())
```

--- 
## 4. How many samples we have with age, antibiotics usage and ctr lpatient

```python
sub = df[df["precise_age_category"].notna()]

age_abx_counts = (
    sub
    .groupby(["precise_age_category", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)

age_abx_counts

sub = df[df["imprecise_age_category"].notna()]

age_abx_counts2 = (
    sub
    .groupby(["imprecise_age_category", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)

age_abx_counts2

```

---

## 5. Age + Sex + antibiotics
```python
age_sex_abx_counts = (
    sub
    .groupby(["precise_age_category", "sex", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)
```
Results:
* age_sex_abx_counts
* Antibiotics_used               No  Yes
* precise_age_category sex              
* Child                female   291   23
*                      male     358   25
* Infant               female   270   81
*                       male     296   52
* Middle-Age Adult     female  1620   44
*                      male    1622   20
* Older Adult          female   699   69
*                      male     925   73
* Oldest Adult         female    82    9
*                      male     100   18
* Teenage              female   543    0
*                      male     579    0
* Toddler              female    40    6
*                      male      52   19
* Unknown              female  2124  196
*                      male    1938  137
* Young adult          female  1235   94
*                      male    1045   90

---

## 5. Age + Sex + Uti + Abx

```python
age_sex_uti_abx_counts = (
    sub
    .groupby(["precise_age_category", "sex", "UTI_history", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)
age_sex_uti_abx_counts
```
  
---

## 6. Age + Sex + BMI + abx 

```pyhton
age_sex_bmi_abx = (
    sub
    .groupby(["precise_age_category", "sex", "BMI_range_new", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)

age_sex_bmi_abx

age_sex_bmi_abx.reset_index().to_csv(
    "age_sex_bmi_antibiotics_counts.tsv",
    sep="\t",
    index=False
)

```
---

## 7. Age + Sex + Gastrointestinal disease + Abx

```python
age_sex_GI_abx = (
    sub
    .groupby(["precise_age_category", "sex", "GI_disease_history", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)

age_sex_GI_abx

age_sex_GI_abx.reset_index().to_csv(
    "age_sex_GI_antibiotics_counts.tsv",
    sep="\t",
    index=False
)
```


---

## 8. Age + Sex + Cancer and adenomas + Abx

```python
Cancers_and_adenomas

age_sex_CA_abx = (
    sub
    .groupby(["precise_age_category", "sex", "Cancers_and_adenomas", "Antibiotics_used"])
    .size()
    .unstack(fill_value=0)
)

age_sex_CA_abx

age_sex_CA_abx.reset_index().to_csv(
    "age_sex_CA_antibiotics_counts.tsv",
    sep="\t",
    index=False
)

```

---

## 9. Preparing the TSE object

```python
# 1. Set acc numbers as the rownames for colData
df = df.set_index("acc", drop=True)
df.index.is_unique

# 2. Replace empty strings with NA
df = df.replace("", pd.NA)

# 3. Drop list-like columns
df = df.drop(columns=["antibiotics_list", "name_of_antibiotic"], errors='ignore')

# 4. Convert categorical columns
categorical_cols = [
    "sex", "age_category", "Antibiotics_used", "UTI_history",
    "HIV_history", "COVID19_history", "CPE_infection_history",
    "Multidrug_resistant_infection", "Antivirals_used"
]
for c in categorical_cols:
    if c in df.columns:
        df[c] = df[c].astype("category")

# 4. Quick checks
print(df.index.is_unique)          # True
print(df.dtypes)                   # Check types
print(df.isna().sum().sort_values(ascending=False).head(10))  # Missing values

5. Drop column subject_disease_status_full

df = df.drop(columns=["subject_disease_status_full"])
```

## 10. Lets modify the antibiotic columns:

```

antibiotic_cols = [f'antibiotic_{i}' for i in range(1, 10)]

all_antibiotics = pd.Series(df[antibiotic_cols].values.ravel()).dropna().unique()

for ab in all_antibiotics:
    df[ab] = 'No'

for col in antibiotic_cols:
    for ab in df[col].dropna().unique():
        df.loc[df[col] == ab, ab] = 'Yes'

df = df.drop(columns=antibiotic_cols, errors='ignore')

df.to_csv("colData_TSE.tsv", sep="\t", index=True)

```

