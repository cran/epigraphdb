% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mappings.R
\name{mappings_gene_to_protein}
\alias{mappings_gene_to_protein}
\title{Return protein uniprot_id from associated genes}
\usage{
mappings_gene_to_protein(
  gene_name_list = NULL,
  gene_id_list = NULL,
  by_gene_id = FALSE,
  mode = c("table", "raw")
)
}
\arguments{
\item{gene_name_list}{List of HGNC symbols of the genes (default)}

\item{gene_id_list}{List of Ensembl gene IDs (when \code{by_gene_id == TRUE})}

\item{by_gene_id}{Search for gene ids (Ensembl gene IDs) instead of gene names (HGNC symbols)}

\item{mode}{If \code{mode = "table"}, returns a data frame
(a \href{https://tibble.tidyverse.org/}{\code{tibble}} as per
\href{https://style.tidyverse.org/}{\code{tidyverse}} convention).
If \code{mode = "raw"}, returns a raw response from EpiGraphDB API
with minimal parsing done by \href{https://httr.r-lib.org/}{\code{httr}}.}
}
\value{
Data from \code{POST /mappings/gene-to-protein}
}
\description{
\href{https://docs.epigraphdb.org/api/api-endpoints/#post-mappingsgene-to-protein}{\code{POST /mappings/gene-to-protein}}
}
\examples{
# By HGNC symbols
\dontrun{
mappings_gene_to_protein(gene_name_list = c("GCH1", "MYOF"))
}

# By Enselbl Ids
\dontrun{
mappings_gene_to_protein(gene_id_list = c("ENSG00000162594", "ENSG00000113302"), by_gene_id = TRUE)
}
}
