# 1. Metadata curation workflow 

This document describes the metadata processing workflow. 

* Input file: Filtered_Metadata.tsv

---

# 1. Download the file

```python
import numpy as np
import pandas as pd
import re


df = pd.read_csv("Filtered_Metadata.tsv", sep="\t")
```
filtered_df.to_csv("Filtered_Metadata.tsv", sep="\t", index=False)
filtered_df = pd.read_csv("Filtered_Metadata.tsv", sep="\t")


## 6. subject_disease_status_full column

### 6.1. inspect the values
```

df["subject_disease_status_full"].value_counts(dropna=False)

#subject_disease_status_full
#NaN                                                            23384
#melanoma                                                         197
#CTR                                                              180
#obesity                                                          166
#iCD                                                              161
#prediabetes                                                       97
#Non-alcoholic fatty liver disease                                 96
#PD                                                                93
#cCD                                                               28
#control                                                           28
#adenocarcinoma                                                    26
#hypertension                                                      23
#adenoma                                                           20
#IC                                                                17
#UC                                                                12
#pre-hypertension                                                  12
#oligohydramnios                                                   11
#Type 2 Diabetes, obesity                                          10
#mild preeclampsia; breech                                         10
#mild preeclampsia                                                 10
#preeclampsia                                                       9
#asthma                                                             4
#Control                                                            4
#Myalgic encephalomyelitis/chronic fatigue syndrome (ME/CFS)        3
#control patient                                                    2
#type 2 diabetes                                                    1
#N                                                                  1

# abbrevations (iCD, cCD, IC and UC) are not clear to interpret so we will search more).
#We may want to use cancers so I'll collect those. Including adenoma (pre-cancer).


```
## 6.2. Cancer related values
```
Cancer_keywords = ["melanoma", "adenocarcinoma"]

df["Cancer"] = pd.NA

def add_label(existing, new):
    if pd.isna(existing):
        return new
    if new in existing:
        return existing
    return f"{existing};{new}"

for keyword in Cancer_keywords:
    mask = df["subject_disease_status_full"].astype(str).str.contains(keyword, case=False, na=False)
    df.loc[mask, "Cancer"] = df.loc[mask, "Cancer"].apply(add_label, new=keyword)

df["Cancer"].value_counts(dropna=False)
```
Cancer
* <NA>              24362
* melanoma            197
* adenocarcinoma       26


```
## 6.2. Abbrevations
```
Diseases_of_interest = ["iCD", "cCD", "IC", "UC", "PD"]

subset = df[df["subject_disease_status_full"].isin(Diseases_of_interest)]

bioprojects = subset["bioproject"].unique()
print("Unique bioproject numbers for iCD, cCD, IC, UC, PD:")
print(bioprojects)

print(bioprojects)
['PRJNA237362' 'PRJDB15956']

## PRJDB15956 is about parkinsos disease (PD)
## PRJDB237362 is about Chrohns disease (CD), Ulcerative colitis (UC)

```python 

IBD_keywords = ["CD", "UC"]

df['IBD'] = df['IBD'].fillna(
    df['subject_disease_status_full'].apply(
        lambda x: "Yes" if x in IBD_keywords else "No"))

df["IBD"].value_counts(dropna=False)
```


## 7. subject_disease_status column`df["subject_disease_status"].value_counts(dropna=False)

```

df["subject_disease_status"].value_counts(dropna=False)

subject_disease_status
#COHORT                                                      10713
#CTR                                                          3651
#NaN                                                          3096
#Crohn's disease                                               936
#Parkinson's disease                                           721
#ulcerative colitis                                            640
#colorectal cancer                                             456
#melanoma                                                      444
#control patient                                               402
#helminthiasis                                                 298
#COVID-19                                                      294
#Clostridium difficile infection                               292
#postmenopausal osteopenia                                     160
#gestational diabetes mellitus                                 133
#type 2 diabetes                                               131
#critically ill patient                                        124
#end-stage renal disease                                       124
#HIV                                                           120
#metabolic syndrome                                            119
#breast cancer                                                 118
#prediabetes                                                    97
#ankylosing spondylitis                                         91
#cholera                                                        89
#metabolic dysfunction-associated steatohepatitis               87
#liver cirrhosis                                                87
#irritable bowel syndrome                                       77
#atopy                                                          73
#postmenopausal osteoporosis                                    72
#nonadvanced colorectal adenoma                                 67
#urinary tract infection                                        67
#autism spectrum disorder                                       66
#konzo                                                          60
#diarrhea                                                       52
#Vogt-Koyanagi-Harada disease                                   51
#polycystic ovary syndrome                                      50
#slow transit constipation                                      48
#chronic kidney disease                                         45
#type 1 diabetes                                                44
#metastatic renal cell carcinoma                                43
#hypertension                                                   34
#impaired glucose tolerance                                     33
#rheumatoid arthritis                                           32
#indeterminate colitis                                          29
#metabolic dysfunction-associated steatotic liver disease       29
#non-small cell lung cancer                                     29
#adenoma                                                        27
#Carbapenemase-producing Enterobacteriaceae infection           27
#chronic obstructive pulmonary disease                          27
#familial adenomatous polyposis                                 22
#Behcet's disease                                               21
#schistosomiasis                                                18
#pre-hypertension                                               12
#stunted growth                                                 12
#short bowel syndrome                                           11
#myalgic encephalomyelitis/chronic fatigue syndrome              3
#advanced colorectal adenoma                                     1
#Name: count, dtype: int64

```
### 7.1. Cancer related values 
```
new_cancer_keywords = [
    "colorectal cancer",
    "melanoma",
    "breast cancer",
    "metastatic renal cell carcinoma",
    "non-small cell lung cancer",
]

def add_label(existing, new):
    if pd.isna(existing):
        return new
    if new in existing:
        return existing
    return f"{existing};{new}"

for keyword in new_cancer_keywords:
    mask = df["subject_disease_status"].astype(str).str.contains(keyword, case=False, na=False)
    df.loc[mask, "Cancer"] = df.loc[mask, "Cancer"].apply(add_label, new=keyword)

df["Cancer"].value_counts(dropna=False)
```

Cancer
* <NA>                                23515
* melanoma                              444
* colorectal cancer                     430
* breast cancer                         118
* metastatic renal cell carcinoma        43
* non-small cell lung cancer             29
* adenocarcinoma;colorectal cancer       26

### 7.2. GI-tract related values 

```
gi_inflam_infect_keywords = [
    "Crohn's disease",
    "ulcerative colitis",
    "indeterminate colitis",]

mask = df['subject_disease_status'].isin(gi_inflam_infect_keywords)
df.loc[mask, 'IBD'] = "Yes"

df["IBD"].value_counts(dropna=False)
```
IBD
* No     23000
* Yes     1663


### 7.3. UTI related values
```
uti_mask = df["subject_disease_status"].astype(str).str.contains("urinary tract infection", case=False, na=False)

df.loc[uti_mask, "UTI_History"] = "Yes"

df["UTI_History"].value_counts(dropna=False)
```

#UTI_history
#No     24355
#Yes      250


### 7.5. Remove original column and save
```
df = df.drop(columns=["subject_disease_status"])
df.to_csv("kesken3.tsv", sep="\t", index=False)
```
## 8 raw_metadata_subject_disease_status_full column
```
df["raw_metadata_subject_disease_status_full"].value_counts(dropna=False)
df = df.drop(columns=["raw_metadata_subject_disease_status_full"])
```
## 9. raw_metadata_host_disease column
```
df["raw_metadata_host_disease"].value_counts(dropna=False)
#raw_metadata_host_disease
#NaN                             24585
#Acute Lymphoblastic Leukemia       14
#Acute Myeloid Leukemia              6
#Name: count, dtype: int64

mask = df["raw_metadata_host_disease"].notna()

df.loc[mask, "Cancer"] = (df.loc[mask, "Cancer"]
    .apply(add_label, new="Leukemia"))

df = df.drop(columns=["raw_metadata_host_disease"])

df.to_csv("kesken3.tsv", sep="\t", index=False)
```
## 10. raw_metadata_diseases column
```
df["raw_metadata_diseases"].value_counts(dropna=False)
df = df.drop(columns=["raw_metadata_diseases"])
```
## 11. raw_metadata_intestinal_disease column

```
df["raw_metadata_intestinal_disease"].value_counts(dropna=False)
df = df.drop(columns=["raw_metadata_intestinal_disease"])
```

## 12. Medications
```python
df["raw_metadata_medications"].value_counts(dropna=False)

df.loc[
    (df['raw_metadata_medications'] == 'antibiotic') &
    (df['Antibiotic_Use'].isna()),
    'Antibiotic_Use'] = 'Yes'

df["Antibiotic_Use"].value_counts(dropna=False)
df = df.drop(columns=["raw_metadata_medications"])
```
Antibiotic_Use
* NaN    21878
* No      1425
* Yes     1302

## medication_with_parents colum
```
df["medication_with_parents"].value_counts(dropna=False)

df = df.drop(columns=["medication_with_parents"])

```
## medication column
```
df["medication"].value_counts(dropna=False)

medication
```
raw_metadata_gestational_age_category
df["raw_metadata_gestational_age_category"].value_counts(dropna=False)

raw_metadata_gestational_age_category
NaN             24573
early term         13
full term          13
late preterm        6



raw_age_sam
df["raw_age_sam"].value_counts(dropna=False)
df = df.drop(columns=["raw_age_sam"])

df = df.drop(columns=["geo_loc_name_country_calc_sra"])
df = df.drop(columns=["raw_metadata_exercise_frequency(n/week)"])
df = df.drop(columns=["raw_metadata_x1_in_the_past_2_weeks_have_you_received_any_of_the_following_medications"])
df = df.drop(columns=["raw_metadata_age_group"])
df = df.drop(columns=["age_group_sam_s_dpl63"])
df = df.drop(columns=["raw_metadata_disease_cause"])
df = df.drop(columns=["raw_metadata_siblings"])

df = df.drop(columns=["study_accession"])

study_accession


print(df.shape)
DataFrame shape: (24605, 39)

df.to_csv("kesken.tsv", sep="\t", index=False)
```

