setwd("~/My Drive/r_packages/bulkPipeline/")
#devtools::load_all()

setwd("~/Desktop/test_bulk/")

df <- pre_process_bulk(counts_filepath = "counts.txt",
                       sample_name_regex = "S[0-9]+",
                       metadata = metadata,
                       edger_min_count = 10,
                       export_normalised_data = TRUE,
                       export_pca = TRUE,
                       export_gene_boxplots = TRUE,
                       boxplot_genes = c("CD3D", "IFNG", "MYC"))

metadata <- filter(metadata, sample %notin% c("s6", "s14"))

contrast_matrix <- makeContrasts(unstim_wt-stim_wt,
                                 unstim_wt-unstim_cpm,
                                 unstim_wt-unstim_c5l2,
                                 stim_wt-stim_cpm,
                                 stim_wt-stim_c5l2,
                                 unstim_cpm-unstim_c5l2,
                                 stim_cpm-stim_c5l2,
                                 levels = metadata$group)

df <- process_bulk(edger_object = df,
                   metadata = metadata,
                   contrast_matrix = contrast_matrix,
                   export_degs = TRUE)

