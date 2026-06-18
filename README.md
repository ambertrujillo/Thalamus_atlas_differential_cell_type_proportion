# Differential Cell Type Proportion Pipeline - Thalamus atlas

## Overview

This pipeline tests for differential cell type proportions (DCTP) across disease conditions using
`propeller` (via the `speckle` package). It consists of two scripts:

1. **`Differential_cell_type_proportion.R`** — transforms proportions, fits the linear model, and outputs results
2. **`Plotting_for_DCTP.R`** — visualizes results as volcano plots and stacked bar plots

---

## Directory Structure

```
Thalamus_atlas_differential_cell_type_proportion/
├── README.md
    ├── scripts
        ├── propeller_analysis.R
        ├── propeller_plotting.R

```

---

## Dependencies

```r
BiocManager::install("speckle")
install.packages(c("limma", "edgeR", "pheatmap", "gt", "anndata",
                   "dplyr", "ggpubr", "ggplot2", "ggrepel"))
```

---

## Script 1: `Differneital_cell_type_proportion.R`

### User Configuration

Within the script, set the following variables before running:

| Variable | Description | Example |
|---|---|---|
| `working_dir` | Path to working directory | `"/path/to/working/directory"` |
| `data_dir` | Path to input annotated `.h5ad` file | `"/path/to/data/directory/adata.h5ad"` |
| `sample` | Column name for donor IDs in `obs` | `"SampleID"` |
| `condition` | Column name for disease/condition in `obs` | `"general_disease"` |
| `cluster` | Column name for cell type annotations in `obs` | `"annotation"` |
| `study` | Column name for batch/cohort blocking variable in `obs` | `"study"` |
| `design` | R formula for design matrix | `~ 0 + condition` |
| `transform` | Proportion transformation | `"asin"` or `"logit"` |
| `contrast_list` | Named list of pairwise contrasts | See below |

### Specifying Contrasts

Contrasts are defined as a named list where each entry is a string specifying the comparison.
Group names must exactly match the factor levels in your condition column after design matrix
cleanup. E.g., A_vs_B = "A - B"; this means that, in the output, +LogFC is higher proportion of a cell type in A and -LogFC means lower proportion of a cell type in A.

```r
contrast_list <- list(
  AD_vs_Control  = "AD - Control",
  AD_vs_FTD      = "AD - FTD",
  FTD_vs_Control = "FTD - Control"
)
```

Each contrast produces one output CSV file named `{contrast_name}_DCTP.csv`.

### Design Formula

The default design `~ 0 + condition` fits a model with one coefficient per condition group,
which is required for pairwise contrasts using `makeContrasts`. To add covariates (e.g., Sex),
modify the formula:

```r
design_formula <- ~ 0 + condition + Sex
```

Note that `study` is handled separately as a **blocking variable** via `duplicateCorrelation`
(random effect) rather than as a fixed covariate in the design matrix, since cohort is treated
as a nuisance variable rather than a variable of interest.

### Output Files

| File | Description |
|---|---|
| `{contrast_name}_DCTP.csv` | One file per contrast with propeller results |
| `Proportions.csv` | Long-format dataframe with raw and transformed proportions per donor per cell type, merged with sample-level metadata |

#### Columns in `{contrast_name}_DCTP.csv`

| Column | Description |
|---|---|
| `logFC` | Log fold change in proportion |
| `AveExpr` | Average transformed proportion across samples |
| `t` | Moderated t-statistic |
| `P.Value` | Raw p-value |
| `adj.P.Val` | BH-adjusted p-value |
| `B` | Log-odds of differential proportion |

#### Columns in `Proportions.csv`

| Column | Description |
|---|---|
| `sample` | Donor ID |
| `clusters` | Cell type label |
| `Proportions` | Raw proportion of that cell type in that donor |
| `TransformedProportions` | Arcsin or logit transformed proportion |
| `condition` | Disease/condition group |
| `study` | Cohort/batch |

---

## Script 2: `Plotting_for_DCTP.R`

Within the script, set the following variables before running:

### User Configuration

| Variable | Description | Example |
|---|---|---|
| `working_dir` | Same working directory as analysis script | `"/data/trujilloae/propeller"` |
| `padj_cutoff` | Adjusted p-value significance threshold for volcano plot | `0.05` |
| `lfc_cutoff` | LogFC threshold for coloring volcano plot points | `0.01` |
| `celltype_colors` | Named vector of colors per cell type | See below |
| `cluster_order` | Order of cell types in stacked bar plot | See below |

### Specifying Colors and Cluster Order

```r
celltype_colors <- c(
  "Neuron"     = "#21409a",
  "Lymphocyte" = "#68151f",
  "Myeloid"    = "#d24a75",
  "Astrocyte"  = "#abd158",
  "OL"         = "#fdb515",
  "OPC"        = "#ed7020",
  "BVC"        = "#74858c",
  "CP"         = "#d5c18e"
)

cluster_order <- c("Neuron", "Lymphocyte", "Myeloid", "Astrocyte", "OL", "OPC", "BVC", "CP")
```

### Plots Produced

**1. Volcano plot per contrast**
Points are colored by direction and significance (up-regulated in orange/red, down-regulated
in blue). Significant cell types (adj.P.Val < `padj_cutoff`) are circled in black and labeled.
Threshold lines are drawn at `padj_cutoff` and `±lfc_cutoff`.

**2. Stacked bar plot of raw proportions per condition**
Raw proportions per condition, colored by cell type. Proportions > 3% are labeled
as percentages in white text.

**3. Stacked bar plot of transformed proportions per condition**
Same layout as above but using arcsin/logit transformed proportions. Proportions > 5%
are labeled. Useful for visualizing the scale on which the statistical model operates.

---

## Statistical Notes

- **Transformation**: arcsin square root (`"asin"`) is the default and is recommended for
  proportion data that includes zeros. Logit transformation (`"logit"`) is available but
  requires all proportions to be strictly between 0 and 1.
- **Blocking**: `study` is included as a random effect via `duplicateCorrelation` to account
  for cohort-level batch effects without losing degrees of freedom to a fixed effect term.
- **Multiple testing**: p-values are adjusted using Benjamini-Hochberg FDR correction within
  each contrast independently via `topTable`.
- **Robust eBayes**: `eBayes(robust = TRUE)` downweights outlier cell types with extreme
  variance, improving stability across cell types with very different proportion distributions.

---

## Running the Pipeline

Interactively in RStudio — set the configuration variables of each script
before sourcing.

---

## Common Issues

| Error | Likely Cause | Fix |
|---|---|---|
| `colnames mismatch` in `stopifnot` | Sample order in `alt_prop` doesn't match `sample_info` | Check for duplicate `SampleID` values or missing samples in metadata |
| `contrast not found` in `makeContrasts` | Group name in contrast string doesn't match design matrix column | Print `colnames(design)` and adjust contrast strings to match exactly |
| `logit` transformation fails | Zeros present in proportion data | Switch to `"asin"` or remove cell types absent from any sample |
| `duplicateCorrelation` warning | Too few samples per study for reliable correlation estimate | Consider collapsing studies or using study as a fixed effect instead |
| Volcano plot points all grey | `significant` column not computed before plotting | Check that `adj.P.Val` and `logFC` columns exist and are numeric in the results CSV |

---

## Citation

Phipson B, Sim A, Hannan E, Martelotto LG, Barnes J, Oshlack A (2022).
*propeller: testing for differences in cell type proportions in single cell data.*
Bioinformatics. doi:10.1093/bioinformatics/btac582