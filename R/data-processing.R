#' Read and process featurecounts data in edgeR
#'
#' This function takes the input from featurecounts,
#' cleans the expression file with better names, creates an edgeR object,
#' filters the data based on min.count, performs TMM-CPM normalisation,
#' and exports boxplots of the data distributions.
#'
#' @param counts_filepath File path to the featurecounts output
#' @param sample_name_regex Unique identifier for the samples in the column names of the featureCounts
#'    output. Uses stringr::str_extract to pull the unique identifiers
#'    e.g. sample_name_regex for a counts.txt file with column names as fcounts_s1_output, fcounts_s2_output,
#'    fcounts_s3_output, fcounts_s4_output would be "s[1-4]". This renders the final column names as s1, s2, s3, s4.
#'    These final column names should match the metadata file.
#' @param metadata Metadata file with sample information in tidy format. Columns should include
#'    "sample" and "group" at minimum. "sample" column must correspond to "sample_name_regex",
#'    and "group" column denotes the comparison groups i.e. stim or unstim. Group is also used to name the
#'    columns in the normalised output.
#' @param edger_min_count Corresponds to the edgeR `min.count` function
#' @param export_normalised_data Should the normalised data be exported?
#' @param export_pca Should the PCA plot be exported?
#' @param pca_dims Figure dimensions for the PCA plot
#' @param export_gene_boxplots Should boxplots be produced for selected genes?
#' @param boxplot_genes Which genes should be plotted?
#' @param boxplot_dims Figure dimensions for the boxplot plot
#' @param exclude_samples Samples to be excluded. Based on metadata$sample column
#'
#' @examples \dontrun{
#'   data <- pre_process_bulk(counts_filepath = "counts.txt",
#'                            sample_name_regex = "S[0-9]+",
#'                            metadata = metadata,
#'                            edger_min_count = 10,
#'                            export_normalised_data = TRUE)
#' }
#'
#' @return edgeR object with norm factors calculated
#'
#' @export


pre_process_bulk <- function(counts_filepath = "featurecounts/counts.txt",
                             sample_name_regex = "s[0-9]+",
                             metadata = NULL,
                             edger_min_count = 10,
                             export_normalised_data = TRUE,
                             export_pca = TRUE,
                             pca_dims = c(4.5, 4.5),
                             export_gene_boxplots = TRUE,
                             boxplot_genes = NULL,
                             boxplot_dims = c(4, 4),
                             exclude_samples = NULL) {

  ##---------- define notin function

  `%notin%` <- Negate(`%in%`)

  ##---------- create output folders

  cat("---------- Creating output_figures directory for output files \n")
  dir.create("r_output/output_figures", recursive = T, showWarnings = F)
  dir.create("r_output/expression_data", showWarnings = F)

  ##---------- read in data

  cat("\n---------- Reading in data \n")

  raw_data <- utils::read.delim(counts_filepath, skip = 1) %>%
    dplyr::select(-c(Chr, Start, End, Strand, Length)) %>%
    dplyr::arrange(Geneid) %>%
    tibble::column_to_rownames("Geneid")

  samples <- colnames(raw_data)

  cat("--- Samples in the counts file are: ", samples, sep = "\n")

  raw_data <- magrittr::set_colnames(raw_data, stringr::str_extract(string = colnames(raw_data),
                                                                    pattern = sample_name_regex)) %>%
    janitor::clean_names() %>%
    dplyr::select(metadata$sample)

  ##---------- excluding samples

  if (!is.null(exclude_samples)) {

    cat("\n---------- Excluding samples")
    cat("--- Samples to be excluded are: ", exclude_samples, sep = "\n")
    raw_data <- dplyr::select(raw_data, -c(dplyr::all_of(exclude_samples)))
    metadata <- dplyr::filter(metadata, sample %notin% exclude_samples)
    warning("Metadata object has been updated internally, but not in the global environment. To update this, run:

            metadata <- filter(metadata, sample %notin% c(", paste0("\"", exclude_samples, "\"", collapse = ", "), "))")

  }

  ##---------- edgeR processing

  cat("\n---------- Creating edgeR object for filtering and normalisation \n")

  groups <- metadata$group
  y <- edgeR::DGEList(counts = raw_data, group = groups)
  keep_genes <- edgeR::filterByExpr(y, min.count = edger_min_count)
  cat("--- Filtering results: \n")
  cat(paste0("- Genes included = ", sum(keep_genes)))
  cat(paste0("\n- Genes excluded = ", length(keep_genes)-sum(keep_genes)), "\n")

  y <- y[keep_genes, , keep.lib.sizes = F]

  y <- edgeR::calcNormFactors(y)
  norm_data <- edgeR::cpm(y) %>% data.frame()

  if (export_normalised_data == TRUE) {
    norm_data <- magrittr::set_colnames(norm_data, make.names(metadata$group, unique = T)) %>%
      data.frame() %>%
      janitor::clean_names()
    utils::write.csv(norm_data, "r_output/expression_data/tmm_cpm_normalised_data.csv")
    cat("--- Normalised data exported to r_output/expression_data/tmm_cpm_normalised_data.csv \n")
  }

  ##---------- plotting some qc

  # raw data

  cat("\n---------- Calculating and plotting expression distributions \n")

  p1 <- raw_data %>%
    magrittr::set_colnames(colnames(norm_data)) %>%
    tidyr::pivot_longer(cols = tidyselect::peek_vars(), names_to = "sample", values_to = "n") %>%
    ggplot2::ggplot(ggplot2::aes(sample, log2(n+1))) +
    ggplot2::geom_boxplot(fill = "gray90", lwd = 0.2) +
    ggplot2::labs(title = "Raw data") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = ggplot2::element_blank())

  # normalised data

  p2 <- norm_data %>%
    tidyr::pivot_longer(cols = tidyselect::peek_vars(), names_to = "sample", values_to = "n") %>%
    ggplot2::ggplot(ggplot2::aes(sample, log2(n+1))) +
    ggplot2::geom_boxplot(fill = "gray90", lwd = 0.2) +
    ggplot2::labs(title = "TMM-CPM normalised data") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = ggplot2::element_blank())

  patchwork::wrap_plots(p1, p2)
  ggplot2::ggsave("r_output/output_figures/boxplot_distributions.png", width = 10, height = 3.5, dpi = 600)
  cat("--- Gene expression distributions exported to r_output/output_figures/boxplot_distributions.png \n")

  ##---------- pca plotting

  if (export_pca == TRUE) {

    cat("\n---------- Calculating and exporting PCA \n")

    ## run pca
    pca_df <- stats::prcomp(t(log2(norm_data + 1)))

    ## get var
    var_data <- pca_df$sdev^2
    var_data <- var_data/sum(var_data)*100
    var_data <- round(var_data, 1)

    ## export pc var plot
    data.frame(var = var_data,
               pc = 1:length(var_data)) %>%
      ggplot2::ggplot(ggplot2::aes(pc, var)) +
      ggplot2::geom_col(fill = "gray90", col = "black", lwd = 0.2) +
      ggplot2::labs(title = "PCA scree plot", x = "PC", y = "% variance explained") +
      ggplot2::scale_x_continuous(breaks = seq(0, length(var_data), 2)) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(fill = NA))

    ggplot2::ggsave(filename = "r_output/output_figures/pca_scree_plot.png", width = 5, height = 3, dpi = 600)
    cat("--- PCA scree plot exported to r_output/output_figures/pca_scree_plot.png \n")

    ## plot pca
    pca_df <- data.frame(pc1 = pca_df$x[, "PC1"],
                         pc2 = pca_df$x[, "PC2"]) %>%
      dplyr::mutate(group = metadata$group)


    ggplot2::ggplot(pca_df, ggplot2::aes(pc1, pc2)) +
      ggplot2::geom_point(shape = 21, size = 4, alpha = 0.8, ggplot2::aes(fill = group)) +
      ggplot2::labs(x = paste0("PC1: ", var_data[1], "%"),
           y = paste0("PC2: ", var_data[2], "%")) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(),
            legend.key = ggplot2::element_blank(),
            aspect.ratio = 1)

    ggplot2::ggsave(filename = "r_output/output_figures/pca_plot.png", width = pca_dims[1], height = pca_dims[2], dpi = 600)
    cat("--- PCA plot exported to r_output/output_figures/pca_plot.png \n")

  }

  ##---------- export boxplots of different genes

  if (export_gene_boxplots == TRUE) {

    cat("\n---------- Plotting and exporting expression data \n")
    cat("--- Gene to use:", paste(boxplot_genes, collapse = ", "), "\n")

    plots <- lapply(boxplot_genes, function(x) {

      norm_data %>%
        tidyr::pivot_longer(cols = tidyselect::peek_vars(), names_to = "sample", values_to = "expression") %>%
        dplyr::mutate(gene = rep(rownames(norm_data), each = ncol(norm_data)),
               group = rep(metadata$group, times = nrow(norm_data))) %>%
        dplyr::filter(gene %in% x) %>%
        ggplot2::ggplot(ggplot2::aes(group, expression)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, ggplot2::aes(fill = group)) +
        ggplot2::geom_jitter(width = 0.3, height = 0, shape = 21, size = 3, ggplot2::aes(fill = group), alpha = 0.8) +
        ggplot2::scale_y_continuous(limits = c(0, NA)) +
        ggplot2::labs(y = "TMM-CPM normalised expression", title = as.character(x)) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
              panel.border = ggplot2::element_rect(fill = NA),
              legend.position = "none",
              axis.title.x = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))

    })

    grDevices::pdf("r_output/output_figures/boxplots.pdf", onefile = TRUE, width = boxplot_dims[1], height = boxplot_dims[2])
    print(plots)
    grDevices::dev.off()

    cat("--- Gene expression plots exported to r_output/output_figures/boxplots.pdf \n \n")

  }

  cat("...done \n \n")

  return(y)

}

#' Process and perform stats testing with edgeR
#'
#' This function takes an edgeR input (typically from `bulkPipeline::pre_process_bulk()`)
#' creates the design matrix, estimates dispersion, and performs DE testing
#'
#' @param edger_object edgeR object with norm factors calculated
#' @param metadata Metadata for comparisons. metadata$group is used for the design matrix
#' @param stats_test What stats test to use. Currently only glmqlf
#' @param contrast_matrix Contrast matrix for comparisons. Created with the `makeContrasts()`
#'    function. See https://github.com/jackbibby1/bulkPipeline for example
#' @param export_degs Should differential expression testing results be exported?
#'
#' @examples \dontrun{
#'   data <- process_bulk(edger_object = y,
#'                        metadata = metadata,
#'                        contrast_matrix = contrast_matrix)
#' }
#'
#' @return edgeR object with norm factors calculated
#'
#' @export



process_bulk <- function(edger_object = NULL,
                         metadata = NULL,
                         stats_test = "glmqlf",
                         contrast_matrix = NULL,
                         export_degs = TRUE) {

  ##---------- prepping data for comparisons
  ## design matrix

  cat("\n---------- Creating the design matrix \n")
  design <- stats::model.matrix(~0 + metadata$group)
  colnames(design) <- unique(metadata$group)
  cat("--- Design matrix to use: \n")
  print(design)

  ## estimate dispersion
  cat("--- Estimating dispersion \n")
  y <- edgeR::estimateDisp(edger_object, design)

  cat("--- Saving BCV plot to r_output/output_figures/bcv_plot.png \n")
  grDevices::png("r_output/output_figures/bcv_plot.png", res = 600, width = 5, height = 5, units = "in")
  edgeR::plotBCV(y)
  grDevices::dev.off()

  if (stats_test == "glmqlf") {

    cat("\n---------- Performing QLF testing \n")

    ## fit model
    cat("--- Fitting model \n")
    fit <- edgeR::glmQLFit(y, design)

    ## decide contrasts
    cat("--- Contrast matrix to use: \n")
    print(contrast_matrix)

    ## do testing

    cat("\n--- Testing with glmQLFTest \n")
    cat("--- Contrasts to make: \n \n")
    diff_ex <- lapply(colnames(contrast_matrix), function(x) {

      print(x)
      edgeR::glmQLFTest(fit, contrast = contrast_matrix[, x]) %>%
        edgeR::topTags(n = Inf) %>%
        data.frame() %>%
        dplyr::mutate(comparison = x,
               direction = dplyr::case_when(FDR < 0.05 & logFC > 0.6 ~ "up",
                                     FDR < 0.05 & logFC < -0.6 ~ "down",
                                     FDR > 0.05 | abs(logFC) < 0.6 ~ "none"))

    })

    names(diff_ex) <- colnames(contrast_matrix)

    if (export_degs == TRUE) {

      cat("\n--- Exporting DEGs to the differential_expression directory \n")

        dir.create("r_output/differential_expression")

      for (i in names(diff_ex)) {

        save_name <- gsub(pattern = " - ", replacement = "_vs_", x = i)
        utils::write.csv(diff_ex[[i]], file = paste0("r_output/differential_expression/", save_name, ".csv"))

      }

    }

  }

  cat("\n ... done")

  return(y)

}












