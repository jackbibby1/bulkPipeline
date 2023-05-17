setwd("~/My Drive/r_packages/bulkPipeline/")
devtools::load_all()

setwd("~/Desktop/test_bulk/")

test <- pre_process_bulk(counts_filepath = "counts.txt",
                         metadata = metadata,
                         boxplot_genes = c("Cpm", "Cd69", "Il2ra"),
                         exclude_samples = c("s6", "s14"))

metadata <- filter(metadata, sample %notin% c("s6", "s14"))

contrast_matrix <- makeContrasts(unstim_wt-stim_wt,
                                 unstim_wt-unstim_cpm,
                                 unstim_wt-unstim_c5l2,
                                 stim_wt-stim_cpm,
                                 stim_wt-stim_c5l2,
                                 unstim_cpm-unstim_c5l2,
                                 stim_cpm-stim_c5l2,
                                 levels = metadata$group)

de_test <- process_bulk(edger_object = test,
                        metadata = metadata,
                        contrast_matrix = contrast_matrix)

