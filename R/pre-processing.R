#' Read and process featurecounts data in edgeR
#'
#' This function takes an input from the featurecounts output,
#' cleans the expression file with better names, creates an edgeR object,
#' filters the data based on min.count, performs TMM-CPM normalisation,
#' and exports boxploits of the data distributions.
#'
#' @param counts_filepath File path to the featurecounts output
#' @param sample_name_regex Regex for sample names to rename the columns
#' @param metadata Metadata file with sample information in tidy format. Columns should include
#'    "sample" and "group" at minimum. "sample" column must correspond to "sample_name_regex",
#'    and "group" column denotes the comparison groups i.e. stim or unstim.
#' @param edger_min_count Corresponds to the edgeR `min.count` function
#' @param export_normalised_data Should the normalised data be exported?
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


pre_process_bulk <- function(counts_filepath = NULL,
                             sample_name_regex = "S[0-9]+",
                             metadata = NULL,
                             edger_min_count = 10,
                             export_normalised_data = TRUE,
                             export_pca = TRUE,
                             pca_dims = c(5, 5),
                             export_gene_boxplots = TRUE,
                             boxplot_genes = NULL,
                             boxplot_dims = c(4, 4),
                             exclude_samples = NULL) {

  ##---------- create output folders

  `%notin%` <- Negate(`%in%`)

  if ("./output_figures" %notin% list.dirs()) {
    cat("---------- Creating output_figures directory for output files \n")
    dir.create("output_figures")
  }

  ##---------- read in data

  cat("\n---------- Reading in data \n")

  raw_data <- read.delim(counts_filepath, skip = 1) %>%
    select(-c(Chr, Start, End, Strand, Length)) %>%
    arrange(Geneid) %>%
    column_to_rownames("Geneid")

  samples <- colnames(raw_data)

  cat("--- Samples in the counts file are: ", samples, sep = "\n")

  raw_data <- set_colnames(raw_data, str_extract(string = colnames(raw_data), pattern = sample_name_regex)) %>%
    janitor::clean_names() %>%
    select(metadata$sample)

  ##---------- excluding samples

  if (!is_null(exclude_samples)) {

    cat("\n---------- Excluding samples")
    cat("--- Samples to be excluded are: ", exclude_samples, sep = "\n")
    raw_data <- select(raw_data, -c(exclude_samples))
    metadata <- filter(metadata, sample %notin% exclude_samples)

  }

  ##---------- edgeR processing

  cat("\n---------- Creating edgeR object for filtering and normalisation \n")

  groups <- metadata$group
  y <- DGEList(counts = raw_data, group = groups)
  keep_genes <- filterByExpr(y, min.count = edger_min_count)
  cat("--- Filtering results: \n")
  cat(paste0("- Genes included = ", sum(keep_genes)))
  cat(paste0("\n- Genes excluded = ", length(keep_genes)-sum(keep_genes)), "\n")

  y <- y[keep_genes, , keep.lib.sizes = F]

  y <- calcNormFactors(y)
  norm_data <- cpm(y) %>% data.frame()

  if (export_normalised_data == TRUE) {
    norm_data <- set_colnames(norm_data, make.names(metadata$group, unique = T)) %>%
      data.frame() %>%
      janitor::clean_names()
    write.csv(norm_data, "tmm_cpm_normalised_data.csv")
  }

  ##---------- plotting some qc

  # raw data

  cat("\n---------- Calculating and plotting expression distributions \n")

  p1 <- raw_data %>%
    set_colnames(colnames(norm_data)) %>%
    pivot_longer(cols = tidyselect::peek_vars(), names_to = "sample", values_to = "n") %>%
    ggplot(aes(sample, log2(n+1))) +
    geom_boxplot(fill = "gray90", lwd = 0.2) +
    labs(title = "Raw data") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank())

  # normalised data

  p2 <- norm_data %>%
    pivot_longer(cols = tidyselect::peek_vars(), names_to = "sample", values_to = "n") %>%
    ggplot(aes(sample, log2(n+1))) +
    geom_boxplot(fill = "gray90", lwd = 0.2) +
    labs(title = "TMM-CPM normalised data") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank())

  patchwork::wrap_plots(p1, p2)
  ggsave("output_figures/boxplot_distributions.png", width = 10, height = 3.5, dpi = 600)
  cat("--- Gene expression distributions exported to output_figures/boxplot_distributions.png \n")

  ##---------- pca plotting

  if (export_pca == TRUE) {

    cat("\n---------- Calculating and exporting PCA \n")

    ## run pca
    pca_df <- prcomp(t(log2(norm_data + 1)))

    ## get var
    var_data <- pca_df$sdev^2
    var_data <- var_data/sum(var_data)*100
    var_data <- round(var_data, 1)

    ## export pc var plot
    data.frame(var = var_data,
               pc = 1:length(var_data)) %>%
      ggplot(aes(pc, var)) +
      geom_col(fill = "gray90", col = "black", lwd = 0.2) +
      labs(title = "PCA scree plot", x = "PC", y = "% variance explained") +
      scale_x_continuous(breaks = seq(0, length(var_data), 2)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA))

    ggsave(filename = "output_figures/pca_scree_plot.png", width = 5, height = 3, dpi = 600)
    cat("--- PCA scree plot exported to output_figures/pca_scree_plot.png \n")

    ## plot pca
    pca_df <- data.frame(pc1 = pca_df$x[, "PC1"],
                         pc2 = pca_df$x[, "PC2"]) %>%
      mutate(group = metadata$group)


    ggplot(pca_df, aes(pc1, pc2)) +
      geom_point(shape = 21, size = 4, alpha = 0.8, aes(fill = group)) +
      labs(x = paste0("PC1: ", var_data[1], "%"),
           y = paste0("PC2: ", var_data[2], "%")) +
      theme(panel.background = element_blank(),
            axis.line = element_line(),
            legend.key = element_blank(),
            aspect.ratio = 1)

    ggsave(filename = "output_figures/pca_plot.png", width = pca_dims[1], height = pca_dims[2], dpi = 600)
    cat("--- PCA plot exported to output_figures/pca_plot.png \n")

  }

  ##---------- export boxplots of different genes

  if (export_gene_boxplots == TRUE) {

    cat("\n---------- Plotting and exporting expression data \n")
    cat("--- Gene to use:", paste(boxplot_genes, collapse = ", "), "\n")

    plots <- lapply(boxplot_genes, function(x) {

      norm_data %>%
        pivot_longer(cols = tidyselect::peek_vars(), names_to = "sample", values_to = "expression") %>%
        mutate(gene = rep(rownames(norm_data), each = ncol(norm_data)),
               group = rep(metadata$group, times = nrow(norm_data))) %>%
        filter(gene %in% x) %>%
        ggplot(aes(group, expression)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.3, aes(fill = group)) +
        geom_jitter(width = 0.3, shape = 21, size = 3, aes(fill = group), alpha = 0.8) +
        labs(y = "TMM-CPM normalised expression", title = as.character(x)) +
        theme(panel.background = element_blank(),
              panel.border = element_rect(fill = NA),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

    })

    pdf("output_figures/boxplots.pdf", onefile = TRUE, width = boxplot_dims[1], height = boxplot_dims[2])
    print(plots)
    dev.off()

    cat("--- Gene expression plots exported to output_figures/boxplots.pdf")

  }

  return(y)

}

#' Process and perform stats testing with edgeR
#'
#' This function takes an edgeR input (typically from `pre_process_bulk()`)
#' creates the design matrix, estimates dispersion, and performs DE testing
#'
#' @param edger_object edgeR object with norm factors calculated
#' @param metadata Metadata for comparisons. metadata$group is used for the design matrix
#' @param stats_test What stats test to use. Currently only glmqlf
#' @param contrast_matrix Contrast matrix for comparisons. Created with the `makeContrasts()`
#'    function
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
  design <- model.matrix(~0 + metadata$group)
  colnames(design) <- unique(metadata$group)
  cat("--- Design matrix to use: \n")
  print(design)

  ## estimate dispersion
  cat("--- Estimating dispersion \n")
  y <- estimateDisp(edger_object, design)

  cat("--- Saving BCV plot to output_figures/bcv_plot.png \n")
  png("output_figures/bcv_plot.png", res = 600, width = 5, height = 5, units = "in")
  plotBCV(y)
  dev.off()

  if (stats_test == "glmqlf") {

    cat("\n---------- Performing QLF testing \n")

    ## fit model
    cat("--- Fitting model \n")
    fit <- glmQLFit(y, design)

    ## decide contrasts
    cat("--- Contrast matrix to use: \n")
    print(contrast_matrix)

    ## do testing

    cat("\n--- Testing with glmQLFTest \n")
    cat("--- Contrasts to make: \n \n")
    diff_ex <- lapply(colnames(contrast_matrix), function(x) {

      print(x)
      glmQLFTest(fit, contrast = contrast_matrix[, x]) %>%
        topTags(n = Inf) %>%
        data.frame() %>%
        mutate(comparison = x,
               direction = case_when(FDR < 0.05 & logFC > 0.6 ~ "up",
                                     FDR < 0.05 & logFC < -0.6 ~ "down",
                                     FDR > 0.05 | abs(logFC) < 0.6 ~ "none"))

    })

    names(diff_ex) <- colnames(contrast_matrix)

    if (export_degs == TRUE) {

      cat("\n--- Exporting DEGs to differential_expression/x... \n")

      if ("./differential_expression" %notin% list.dirs()) {
        dir.create("differential_expression")
      }

      for (i in names(diff_ex)) {

        save_name <- gsub(pattern = " - ", replacement = "_vs_", x = i)
        write.csv(diff_ex[[i]], file = paste0("differential_expression/", save_name, ".csv"))

      }

    }

  }

}












