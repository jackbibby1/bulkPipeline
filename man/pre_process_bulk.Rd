% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre-processing.R
\name{pre_process_bulk}
\alias{pre_process_bulk}
\title{Read and process featurecounts data in edgeR}
\usage{
pre_process_bulk(
  counts_filepath = NULL,
  sample_name_regex = "S[0-9]+",
  metadata = NULL,
  edger_min_count = 10,
  export_normalised_data = TRUE,
  export_pca = TRUE,
  pca_dims = c(4.5, 4.5),
  export_gene_boxplots = TRUE,
  boxplot_genes = NULL,
  boxplot_dims = c(4, 4),
  exclude_samples = NULL
)
}
\arguments{
\item{counts_filepath}{File path to the featurecounts output}

\item{sample_name_regex}{Unique identifier for the samples in the column names of the featureCounts
output. Uses stringr::str_extract to pull the unique identifiers
e.g. sample_name_regex for a counts.txt file with column names as fcounts_s1_output, fcounts_s2_output,
fcounts_s3_output, fcounts_s4_output would be "s[1-4]". This renders the final column names as s1, s2, s3, s4.
These final column names should match the metadata file.}

\item{metadata}{Metadata file with sample information in tidy format. Columns should include
"sample" and "group" at minimum. "sample" column must correspond to "sample_name_regex",
and "group" column denotes the comparison groups i.e. stim or unstim. Group is also used to name the
columns in the normalised output.}

\item{edger_min_count}{Corresponds to the edgeR `min.count` function}

\item{export_normalised_data}{Should the normalised data be exported?}

\item{export_pca}{Should the PCA plot be exported?}

\item{pca_dims}{Figure dimensions for the PCA plot}

\item{export_gene_boxplots}{Should boxplots be produced for selected genes?}

\item{boxplot_genes}{Which genes should be plotted?}

\item{boxplot_dims}{Figure dimensions for the boxplot plot}

\item{exclude_samples}{Samples to be excluded. Based on metadata$sample column}
}
\value{
edgeR object with norm factors calculated
}
\description{
This function takes an input from the featurecounts output,
cleans the expression file with better names, creates an edgeR object,
filters the data based on min.count, performs TMM-CPM normalisation,
and exports boxploits of the data distributions.
}
\examples{
\dontrun{
  data <- pre_process_bulk(counts_filepath = "counts.txt",
                           sample_name_regex = "S[0-9]+",
                           metadata = metadata,
                           edger_min_count = 10,
                           export_normalised_data = TRUE)
}

}
