###################################################################
##### Running Differential cell type proportions -- Propeller #####
###################################################################

working_dir=""
data_dir=""

setwd(working_dir)

#BiocManager::install("speckle")

library(speckle)
library(limma)
library(edgeR)
library(anndata)
library(dplyr)

##### Load data - and make dataframe #####
data <- read_h5ad("adata.h5ad")


data_df <- data.frame(
  "sample" = factor(as.character(data$obs$SampleID)), # Donor IDs
  "condition" = factor(as.character(data$obs$general_disease)), # Condition to be compared
  "cluster" = factor(as.character(data$obs$annotation)), # Clusters you want cell type proportions for
  "study" = factor(as.character(data$obs$study)) # Any other covariates to be included
)

table(data_df$cluster)

##### Transform proportions #####
alt_prop = getTransformedProps(clusters = data_df$cluster, sample = data_df$sample, transform = "asin") # This is an arcsin transformation but logit transformation is also available

##### Set up design matrix #####
sample_info = data_df %>%
  group_by(sample) %>%
  summarise(
    condition = unique(condition),
    study = unique(study)
  ) %>%
  arrange(sample)

# For differential cell type proportion the column names of transformed proportion dataframe need to match the order of the sample name in the metadata
stopifnot(all(colnames(alt_prop$TransformedProps) == sample_info$sample))

design = model.matrix(~ 0 + condition, data=sample_info)
colnames(design) <- gsub("condition", "", colnames(design))  # Clean up names

### This is to create a blocking variable (random effect) in the model
dupcor = duplicateCorrelation(alt_prop$TransformedProps, design=design, block=sample_info$study)

### Fit linear model with random effect
fit <- lmFit(
  alt_prop$TransformedProps, 
  design,
  block = sample_info$study,
  correlation = dupcor$consensus.correlation
)

### Define contrasts for all pairwise comparisons
contrasts <- makeContrasts(
  AD_vs_Control = AD - Control,
  AD_vs_FTD = AD - FTD,
  Control_vs_FTD = FTD - Control,
  levels = design
)

### Fit contrasts
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, robust = TRUE)

### Extract results for each comparison
results_ADvControl <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
results_ADvFTD <- topTable(fit2, coef = 2, number = Inf, sort.by = "P")
results_FTDvControl <- topTable(fit2, coef = 3, number = Inf, sort.by = "P")

### Write results
write.csv(results_ADvControl, file="AD_vs_Control_DCTP.csv")
write.csv(results_ADvFTD, file="AD_vs_FTD_DCTP.csv")
write.csv(results_FTDvControl, file="FTD_vs_Control_DCTP.csv")


### Calculate and save proportions 
props = data.frame(alt_prop$Proportions)
colnames(props)[3] = "Proportions"

transformed_props = data.frame(alt_prop$TransformedProps)
colnames(transformed_props)[3] = "TransformedProportions"

proportions = left_join(props, transformed_props, by=c("sample", "clusters"))
proportions = left_join(proportions, sample_info, by="sample")

write.csv(proportions, file="Proportions.csv", row.names = FALSE)














