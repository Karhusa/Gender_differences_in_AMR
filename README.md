# Master's thesis: Gender differences in ARG load

Saana Karhula 29.10.2025 -

## 0. Data

Includes large datafiles used for combining metadata

* SRA:SRA_metadata_with_biosample.tsv
* Metalog: human_all_wide_2025-10-19.tsv.gz
* merged_all.tsv
* TreeSummarisedExperiment(TSE): TSE.rds
  * Table was too large for github
  
## 1. Table fixing

00_

Matching Metalog to TSE Object, and SRA Metadata by sample identifiers (accession numbers)

Overview:
This workflow aims to determine whether human samples listed in Metalog are represented in an existing TSE object, and to annotate SRA metadata accordingly.
Because Metalog sample lists are very large, Unix tools are used for preprocessing, followed by matching in R.

Result: Matches can't be foumd by using accession numbers.

01_

Processing SRA Metadata: Extracting Gender and Sample Type

Overview:
This workflow extracts human-related metadata from a large SRA attributes file, focusing on the following keywords:

* Sex / Gender (female, male)
* Disease
* Age
* BMI
* Antibiotics
* Sample type (feces / gut)
* The raw metadata contains a JSON-like column (jattr) that must be unpacked before filtering.

02_

Linking SRA Metadata with Metalog Human Metadata via BioSample IDs

Overview:
This workflow links SRA run-level metadata to SPIRE human metadata using shared BioSample IDs (SAMN, SAMD, SAMEA).
The goal is to produce a clean, merged metadata table enriched with demographic, clinical, and antibiotic-related information.

03_ Metadata Cleaning and Harmonization Workflow

Input:
cleaned_merged_final.tsv
(merged SPIRE + SRA metadata, 24,605 rows)

Goal:
Create a clean, analysis-ready metadata table with harmonized BMI, age, disease status, infections, and antibiotic usage.

04


05

---

## 2. TSE




