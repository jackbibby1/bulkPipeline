# bulkPipeline

- Pipeline for basic bulk RNA-seq analysis
- Follows a general featureCounts > edgeR pipeline

#### Installation:

```
devtools::install_github("jackbibby1/bulkPipeline")
```

#### Main functions:

- `pre_process_bulk()` covers: 
  - Reading in featureCounts data
  - Cleaning up the columns
  - Order the columns based on the metadata file
  - Import counts and groups into edgeR
  - Filter genes based on edgeR min.count
  - Normalisation and export of normalised data
  - Exporting some basic QC on count distribution
  - PCA of normalised data
  - Exporting boxplots of selected genes across samples
- `process_bulk()` covers:
  - Creating a design matrix based on metadata groups
  - Estimating dispersion
  - Exporting BCV plot from edgeR
  - Model fitting and testing through glmQLFit > glmQLFTest
  - Exporting DEGs for all comparisons in the design

#### Example:

```ruby

# set wdir and generate metadata ------------------------------------------
setwd("path-to/wdir")
metadata <- data.frame(sample = paste0("sample", 1:12), 
                       stimulation = rep(c("unstim", "stim"), each = 6),
                       disease_status = rep(c("healthy", "disease"), each = 3, times = 2)) %>%
  mutate(group = paste(stimulation, disease_status, sep = "_"))

# pre-process data --------------------------------------------------------
# returns an edgeR object containing cleaned up names, groups, and norm factors calculated
# also generates an "output_figures" directory to store all the plots (QC, PCA, boxplots etc.)

df <- pre_process_bulk(counts_filepath = "counts.txt",
                       sample_name_regex = "S[0-9]+",
                       metadata = metadata,
                       edger_min_count = 10,
                       export_normalised_data = TRUE,
                       export_pca = TRUE,
                       export_gene_boxplots = TRUE,
                       boxplot_genes = c("CD3D", "IFNG", "MYC"))

# create a contrast matrix for comparisons --------------------------------

contrast_matrix <- makeContrasts(unstim_healthy-stim_healthy,
                                 unstim_disease-stim_disease,
                                 unstim_healthy-unstim_disease,
                                 stim_healthy-stim_disease,
                                 levels = metadata$group)

# downstream processing ---------------------------------------------------
# returns an edgeR object with norm factors, dispersions etc. calculated
# also generates a "differential_expression" directory with all the results for each comparison

df <- process_bulk(edger_object = df,
                   metadata = metadata,
                   contrast_matrix = contrast_matrix, 
                   export_degs = TRUE)
                    
```
