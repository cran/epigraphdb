% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pqtl.R
\name{pqtl}
\alias{pqtl}
\title{Return information related to the pQTL analysis}
\usage{
pqtl(
  query,
  rtype = c("mrres", "simple", "sglmr", "inst", "sense"),
  pvalue = 0.05,
  searchflag = c("traits", "proteins"),
  mode = c("table", "raw")
)
}
\arguments{
\item{query}{(Required) A protein coding gene name or a trait name,
e.g. "ADAM19" or "Inflammatory bowel disease"
which cannot be \code{NULL}.}

\item{rtype}{(Optional) A type of data to be extracted, which can be one of these options:
\enumerate{
\item \code{simple}: Basic summary
\item \code{mrres}: MR results (DEFAULT)
\item \code{sglmr}: Single SNP MR results
\item \code{inst}: SNP information
\item \code{sense}: Sensitivity analysis
\strong{NOTE}: \code{mrres} is a DEFAULT option.
}}

\item{pvalue}{(Optional) A pvalue threshold for MR results with the DEFAULT set to 0.05.
\strong{NOTE}: this threshold applies to any \code{rtype} chosen.}

\item{searchflag}{(Required) A flag to indicate whether you are searching for proteins or
traits which cannot be \code{NULL}.
If \code{query} is a protein name, then this flag should be "proteins";
if \code{query} is a trait, this flag should be "traits".
\strong{NOTE}: if the wrong flag is chosen for \code{query}, there will be no result
returned.}

\item{mode}{If \code{mode = "table"}, returns a data frame
(a \href{https://tibble.tidyverse.org/}{\code{tibble}} as per
\href{https://style.tidyverse.org/}{\code{tidyverse}} convention).
If \code{mode = "raw"}, returns a raw response from EpiGraphDB API
with minimal parsing done by \href{https://httr.r-lib.org/}{\code{httr}}.}
}
\value{
Data from \verb{GET /pqtl/}
}
\description{
\href{https://docs.epigraphdb.org/api/api-endpoints/#get-pqtl}{\verb{GET /pqtl/}}
}
\examples{
# Returns a data frame of MR results, while searching for proteins
\dontrun{
pqtl(query = "ADAM19", searchflag = "proteins")
}

# Returns a data frame with SNP information, while searching for traits
\dontrun{
pqtl(
  query = "Inflammatory bowel disease",
  rtype = "inst",
  searchflag = "traits"
)
}

# Change a pvalue threshold (the default is 0.05)
\dontrun{
pqtl(
  query = "Inflammatory bowel disease",
  rtype = "inst",
  pvalue = 1.0,
  searchflag = "traits"
)
}

# Returns raw response if mode="raw"
\dontrun{
pqtl(
  query = "ADAM19", searchflag = "proteins",
  mode = "raw"
) \%>\% str()
}
}
