# Psoriasis Cosmeceutical – Bioinformatics Validation

**Goal**: Strengthen the undergraduate cosmeceutical formulation (CBD + lutein + tulsi + aloe + shea + cetyl alcohol) with **public transcriptomics**, **mechanism mapping**, **connectivity mapping**, and **ADMET/permeation** to produce a **publishable, reproducible** computational study.

## Primary Hypothesis
The actives counter the **psoriatic skin transcriptional program** (esp. IL-17/TNF/NF-κB, keratinocyte hyperproliferation/barrier disruption) and are **computationally consistent** with topical delivery.

## Endpoints
- **Primary**: Negative connectivity between psoriatic disease signature and drug/mechanism signatures; significant **pathway opposition**.
- **Secondary**: ADMET/skin-permeation plausibility; docking sanity checks; (optional) single-cell localization.

## Repro Quickstart

### 1) Create conda env
```bash
conda env create -f env/environment.yml
conda activate psoria-bioinf
```

### 2) Download public datasets
Edit `data/raw/DATASETS.md` with selected GEO accessions, then run:
```bash
python src/utils/download_data.py --accessions GSEXXXXX GSEYYYYY
```

### 3) Run analyses (examples)
```bash
# Bulk DE + GSEA (R)
Rscript src/analysis/01_DE_and_GSEA.R --data data/raw --out results

# Target mapping + PPI (Python)
python src/analysis/02_target_mapping_ppi.py --de results/tables/de_top.csv --out results

# Connectivity mapping (Python)
python src/analysis/03_connectivity_mapping.py --up results/tables/up_genes.txt --down results/tables/down_genes.txt --out results

# ADMET + Skin Permeation (Python + RDKit)
python src/analysis/04_admet_skin.py --compounds data/raw/compounds.csv --out results
```

### 4) Manuscript assembly
See `manuscript/outline.md` for figure callouts (F1–F5) and narrative skeleton.

## Repo Map
```
psoriasis-cosmeceutical-bioinf/
├─ README.md
├─ .gitignore
├─ env/environment.yml
├─ data/
│  ├─ raw/              # GEO accessions + downloads
│  └─ processed/        # intermediate (gitignored)
├─ notebooks/           # optional exploratory
├─ src/
│  ├─ utils/            # download/helpers
│  └─ analysis/         # analysis scripts
├─ results/
│  ├─ figures/
│  └─ tables/
└─ manuscript/
   └─ outline.md
```

## Progress Log
Track work in `PROGRESS_LOG.md` (append entries).

---

**Provenance**: This repo integrates and extends an undergraduate thesis formulation via public bioinformatics validation.