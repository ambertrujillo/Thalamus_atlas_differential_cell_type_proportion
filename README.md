# Differential Cell Type Proportion - Thalamus atlas
---

## Dependencies

```r
BiocManager::install("speckle")
install.packages(c("limma", "edgeR", "pheatmap", "gt", "anndata", "dplyr", "ggpubr", "ggplot2", "ggrepel"))
```

---

## Script 1: `propeller_analysis.R`

### User Configuration

At the top of the script, set the following variables before running:

| Variable | Description | Example |
|---|---|---|
| `working_dir` | Path to working directory | `"/data/trujilloae/propeller"` |
| `data_path` | Path to input `.h5ad` file | `"adata.h5ad"` |
| `sample_col` | Column name for donor IDs in `obs` | `"SampleID"` |
| `condition_col` | Column name for disease/condition in `obs` | `"general_disease"` |
| `cluster_col` | Column name for cell type annotations in `obs` | `"annotation"` |
| `study_col` | Column name for batch/cohort blocking variable in `obs` | `"study"` |
| `design_formula` | R formula for design matrix | `~ 0 + condition` |
| `transform_method` | Proportion transformation | `"asin"` or `"logit"` |
| `contrast_list` | Named list of pairwise contrasts | See below |

### Specifying Contrasts

Contrasts are defined as a named list where each entry is a string specifying the comparison.
Group names must exactly match the factor levels in your condition column after design matrix
cleanup (spaces and special characters are removed automatically).

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

### Output

One CSV per contrast saved to `working_dir`, with the following columns:

| Column | Description |
|---|---|
| `logFC` | Log fold change in proportion |
| `AveExpr` | Average expression (mean transformed proportion) |
| `t` | Moderated t-statistic |
| `P.Value` | Raw p-value |
| `adj.P.Val` | BH-adjusted p-value |
| `B` | Log-odds of differential proportion |

---

## Script 2: `propeller_plotting.R`

### User Configuration

| Variable | Description | Example |
|---|---|---|
| `working_dir` | Same working directory as analysis script | `"/data/trujilloae/propeller"` |
| `results_dir` | Path to results CSVs | `"results/"` |
| `contrast_list` | Same contrast names used in analysis script | See above |
| `padj_cutoff` | Adjusted p-value significance threshold | `0.05` |
| `condition_colors` | Named vector of colors per condition | See below |

### Specifying Colors

```r
condition_colors <- c(
  "Control" = "#3B6D11",
  "AD"      = "#A32D2D",
  "FTD"     = "#534AB7"
)
```

### Plots Produced

**1. Bar plot of cell type proportions per condition**
Shows mean proportion ± SEM per cell type, colored by condition. Significant cell types
(adj.P.Val < `padj_cutoff`) are annotated with asterisks.

**2. Heatmap of transformed proportions**
Rows are cell types, columns are donors, colored by transformed proportion value.
Donors are annotated by condition and study.

**3. Volcano plot per contrast**
x-axis: logFC in proportion; y-axis: -log10(adj.P.Val). Significant cell types are
labeled and colored.

---

## Statistical Notes

- **Transformation**: arcsin square root (`"asin"`) is the default and is recommended for
  proportion data. Logit transformation (`"logit"`) is available as an alternative but
  requires all proportions to be strictly between 0 and 1 (no zeros).
- **Blocking**: `study` is included as a random effect via `duplicateCorrelation` to account
  for cohort-level batch effects without losing degrees of freedom to a fixed effect term.
- **Multiple testing**: p-values are adjusted using Benjamini-Hochberg FDR correction within
  each contrast independently.
- **Robust eBayes**: `eBayes(robust = TRUE)` is used to downweight outlier cell types with
  extreme variance, improving stability across cell types with very different proportion
  distributions.

---

## Running the Pipeline

```bash
# Step 1: Run analysis
Rscript propeller_analysis.R

# Step 2: Run plotting
Rscript propeller_plotting.R
```

Or interactively in RStudio — set the configuration variables at the top of each script
before sourcing.

---

## Common Issues

| Error | Likely Cause | Fix |
|---|---|---|
| `colnames mismatch` in `stopifnot` | Sample order in `alt_prop` doesn't match `sample_info` | Check for duplicate `SampleID` values or missing samples in metadata |
| `contrast not found` in `makeContrasts` | Group name in contrast string doesn't match design matrix column | Print `colnames(design)` and adjust contrast strings to match exactly |
| `logit` transformation fails | Zeros present in proportion data | Switch to `"asin"` or remove cell types with zero counts in any sample |
| `duplicateCorrelation` warning | Too few samples per study for reliable correlation estimate | Consider collapsing studies or using study as a fixed effect instead |

---

## Citation

Phipson B, Sim A, Hannan E, Martelotto LG, Barnes J, Oshlack A (2022).
*propeller: testing for differences in cell type proportions in single cell data.*
Bioinformatics. doi:10.1093/bioinformatics/btac582