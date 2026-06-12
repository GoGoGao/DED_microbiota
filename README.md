# Ocular Surface Microbiome — ECST Analysis Code

This repository provides the reproduction code for the companion paper, covering **ECST (Eye Surface Microbiome Community State Types)** construction, cross-cohort filtering and visualization, co-occurrence network comparison, clinical association analysis, and PAM benchmark validation.

Each plotting script simultaneously outputs the corresponding **TSV source data**; composite / multi-panel figures save individual images and TSV files for each sub-panel.

---

## Directory Structure

```
.
├── input/                              # Shared input data
│   ├── species.relative.abundance.tsv   # Species relative abundance (main analysis, Bradyrhizobium excluded)
│   └── species.relative.abundance.all.tsv  # Full species abundance table (all species included)
│
├── Figure2-species_abundance_distribution/   # Figure 2: Cross-cohort species abundance distribution & filtering
│   ├── cross_cohort_filter_viz.py
│   ├── work.sh
│   └── results/
│
├── Figure3-5AB-constructECST/              # Figure 3–5A/B: ECST construction & networks
│   ├── ECST_pipeline.R                   # Main typing pipeline (M01–M08)
│   ├── ECST_network_compare.R              # Figure 4C/D co-occurrence networks
│   ├── group.tsv                           # Sample ECST assignment results (sample_id × ECST)
│   ├── abundance_corrected.tsv -> ...      # MMUPHin-corrected abundance (symlink)
│   ├── work.sh
│   ├── results/                            # Main ECST pipeline outputs
│   └── result_network/                     # Network analysis outputs
│
├── Figure5CtoJ-ClinicalFactors/            # Figure 5C–J: Clinical factor associations
│   ├── ecst_clinical_figure.py
│   ├── clinical.tsv                        # Clinical data + ECST assignments
│   ├── work.sh
│   └── output/
│
└── Figure6-benchmarkPAM/                   # Figure 6: PAM vs DMM benchmark comparison
    ├── pam.r
    ├── work.sh
    └── result/
```

---

## Environment Requirements

A Conda environment is recommended:

```bash
conda activate eyemicrobiota_env
```

### Python Dependencies

```bash
pip install pandas numpy matplotlib scipy kneed openpyxl
```

| Module | Script |
|--------|--------|
| Figure 2 | `cross_cohort_filter_viz.py` |
| Figure 5C–J | `ecst_clinical_figure.py` |

### R Dependencies

| Package | Purpose |
|---------|---------|
| MMUPHin | Batch effect correction |
| DirichletMultinomial | DMM typing |
| vegan, cluster, ape | Diversity / clustering / ordination |
| tidyverse, patchwork, ggpubr, ggrepel | Plotting |
| ComplexHeatmap, circlize, gridExtra | Heatmaps |
| slingshot, SingleCellExperiment | Trajectory analysis (optional) |
| WGCNA, igraph, ggraph | Co-occurrence networks |
| optparse | Command-line argument parsing |

Installation example (Bioconductor + CRAN):

```bash
# Run inside R
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("MMUPHin", "DirichletMultinomial", "ComplexHeatmap",
                        "slingshot", "SingleCellExperiment"))
install.packages(c("vegan", "cluster", "tidyverse", "patchwork", "ggpubr",
                 "ggrepel", "circlize", "gridExtra", "WGCNA", "igraph",
                 "ggraph", "ggsci", "scales", "optparse"))
```

---

## Input Data Description

### 1. Species Relative Abundance Tables (`input/`)

| File | Description |
|------|-------------|
| `species.relative.abundance.tsv` | Main analysis abundance table — *Bradyrhizobium* **excluded** |
| `species.relative.abundance.all.tsv` | Complete abundance table — **all** species included |

Format requirements:

- TSV; first column is species name (Kraken2 format `k__Bacteria|s__Species_name` supported); remaining columns are sample IDs.
- Values are relative abundances (0–1 or 0–100; the scripts auto-detect the scale).

### 2. Sample Metadata (to be provided by the user)

Required metadata columns for each module:

| Column | Figure 2 | ECST Pipeline | PAM Benchmark |
|--------|:--------:|:-------------:|:-------------:|
| `Sample` | ✓ | ✓ | ✓ |
| `Dataset` (cohort / dataset) | ✓ | ✓ | ✓ |
| `Group` (disease / control) | optional | ✓ | ✓ |
| `Country` | optional | optional | ✓ |
| `SamplingSite` | optional | optional | ✓ |

> **Note:** `Figure3-5AB-constructECST/group.tsv` in this repository contains ECST assignment results (`sample_id` + `ECST`) and **cannot** be used as a substitute for the metadata above. Before re-running the analysis, save the complete metadata as `input/metadata.tsv` (or update the `--metadata` paths in each `work.sh` as needed).

### 3. Clinical Data (`Figure5CtoJ-ClinicalFactors/clinical.tsv`)

This file integrates ECST assignments with clinical indicators, including columns: `sample`, `ECST`, `Group`, `MGD_grade`, `Age`, `OSDI`, `Schirmer`, `TMH`, `MGL`, `TBUT`, etc. For reproducing Figure 5C–J only, this file can be used directly without additional input.

---

## Recommended Execution Order

```
Prepare input data
    │
    ├─► Figure 2 (cross-cohort filtering visualization, independent)
    │
    └─► Figure 3–5A/B main pipeline (ECST_pipeline.R)
            │
            ├─► Generates group.tsv + abundance_corrected.tsv
            ├─► ECST_network_compare.R (Figure 4C/D)
            │
            ├─► Figure 5C–J (clinical associations, uses clinical.tsv)
            │
            └─► Figure 6 (PAM benchmark validation)
```

---

## Module Usage

### Figure 2 — Cross-Cohort Species Abundance Distribution (Figure 2A–E)

Cross-cohort consistency filtering and Kneedle elbow detection, visualizing the abundance distribution of retained vs. filtered species.

```bash
cd Figure2-species_abundance_distribution
bash work.sh
```

Or run manually:

```bash
python cross_cohort_filter_viz.py \
    --abundance ../input/species.relative.abundance.tsv \
    --metadata  ../input/metadata.tsv \
    --outdir    results/ \
    --min-cohorts 2 \
    --top-n 10 \
    --kneedle-range 5.0 \
    --dpi 300
```

**Outputs** (`results/`):

| File | Corresponding Sub-figure |
|------|--------------------------|
| `fig_A_abundance_distribution.svg/.tsv` | Figure 2A: Per-cohort species abundance distribution |
| `fig_B_kneedle_elbow.svg/.tsv` | Figure 2B: Kneedle elbow curve |
| `fig_C_retained_stacked_bar.svg/.tsv` | Figure 2C: Retained core species stacked bar chart |
| `fig_D_filtered_stacked_bar.svg/.tsv` | Figure 2D: Filtered cohort-specific species |
| `fig_E_retained_profile.svg/.tsv` | Figure 2E: Retained species abundance profiles |
| `S_tables.xlsx` | Excel summary of all TSV supplementary tables |

---

### Figure 3–5A/B — ECST Construction (Figure 3A–G, 4A–B, 5A–B)

#### Step 1: Main Typing Pipeline

```bash
cd Figure3-5AB-constructECST

Rscript ECST_pipeline.R \
    --abundance ../input/species.relative.abundance.tsv \
    --metadata  ../input/metadata.tsv \
    --outdir    results \
    --manual_k  4
```

Main steps: species name parsing & deduplication → cross-cohort pre-filtering → MMUPHin batch correction → DMM typing (k=4) → diversity / ordination / trajectory / validation.

**Optional Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mean_abd_thr` | 0.0054 | Within-cohort mean abundance threshold (0.54%) |
| `--min_cohorts` | 2 | Minimum number of cohorts a species must pass in |
| `--manual_k` | 4 | Number of ECST types; set to 0 for automatic selection |
| `--kmax` | 8 | Maximum k scanned by DMM |
| `--exclude` | — | Cohorts to exclude, comma-separated |
| `--seed` | 2025 | Random seed |

**Output Directories** (`results/`):

| Subdirectory | Contents |
|--------------|----------|
| `M01_MMUPHin/` | Pre/post-correction PCoA, PERMANOVA, `abundance_corrected.tsv` |
| `M02_DMM/` | DMM model selection, Silhouette, `ecst_assignment.tsv` |
| `M03_Composition/` | Composition heatmap, ECST naming |
| `M04_Ordination/` | PCoA coordinates, Slingshot pseudotime |
| `M05_Diversity_ECST/` | Diversity comparisons (Figure 3A–D, 4A) |
| `M07_Disease/` | Disease proportions (Figure 5A–B) |
| `M08_Validation/` | LOO/LOCO validation, core species consistency (Figure 3E–G) |

Each subdirectory contains the corresponding **PDF/PNG figures** and **TSV source data**.

#### Step 2: Extract ECST Assignments and Prepare Network Analysis

```bash
# Extract sample_id + ECST from DMM results
cut -f1,2 results/M02_DMM/ecst_assignment.tsv > group.tsv

# Symlink the corrected abundance table
ln -sf results/M01_MMUPHin/abundance_corrected.tsv .

# Run co-occurrence network comparison (Figure 4C/D)
Rscript ECST_network_compare.R
```

Or run all steps at once:

```bash
bash work.sh
```

**Network Outputs** (`result_network/`):

| Path | Description |
|------|-------------|
| `plots/Figure4C_*.pdf/.png/.tsv` | Per-ECST co-occurrence network panels |
| `plots/Figure4D_*.pdf/.png/.tsv` | Network topology metric bar charts |
| `summary/network_topology_comparison.tsv` | Topology statistics summary |
| `ECST_*/edges.tsv`, `node_properties.tsv` | Network edges and node properties per ECST |

---

### Figure 5C–J — Clinical Factor Associations

Generates plots associating ECST with clinical indicators (OSDI, MGD, MGL, TBUT, TMH, Age, etc.) based on `clinical.tsv`.

```bash
cd Figure5CtoJ-ClinicalFactors
bash work.sh
```

Or:

```bash
python ecst_clinical_figure.py --input clinical.tsv --outdir ./output
```

**Outputs** (`output/`):

| Figure | TSV Supplementary Table |
|--------|------------------------|
| `Figure5C.pdf/.jpg/.svg` | `STable_Figure5C.tsv` |
| `Figure5D` – `Figure5J` | `STable_Figure5D.tsv` – `STable_Figure5J.tsv` |
| — | `STable_Figure5.xlsx` (all supplementary tables combined) |

---

### Figure 6 — PAM Benchmark Validation (Figure 6A–D)

Compares PAM clustering (k=3) with DMM typing (k=4) to evaluate clustering stability and disease associations.

```bash
cd Figure6-benchmarkPAM
bash work.sh
```

Or:

```bash
Rscript pam.r \
    --abundance ../input/species.relative.abundance.tsv \
    --metadata  ../input/metadata.tsv \
    --outdir    result \
    --manual_k  4 \
    --pam_k     3
```

**Outputs** (`result/`):

| Subdirectory / File | Corresponding Figure |
|---------------------|----------------------|
| `M03_PAM_Validation/Figure6A_Silhouette.*` + `.tsv` | Figure 6A: Silhouette comparison |
| `M03_PAM_Validation/Figure6B_PCoA_PAM.*` + `.tsv` | Figure 6B: PAM PCoA |
| `M05_PAM_Validation/Figure6C_LOCO_ARI.*` + `.tsv` | Figure 6C: LOCO ARI |
| `M04_PAM_Disease/Figure6D_Disease_Proportion.*` + `.tsv` | Figure 6D: Disease proportions |

---

