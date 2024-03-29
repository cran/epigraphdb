% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.R
\name{meta_nodes_list_node}
\alias{meta_nodes_list_node}
\title{List nodes under a meta node}
\usage{
meta_nodes_list_node(
  meta_node,
  full_data = TRUE,
  limit = 10,
  offset = 0,
  mode = c("table", "raw")
)
}
\arguments{
\item{meta_node}{Name of a meta node (e.g. Gwas). Use \code{meta_nodes_list} to get the full list of meta nodes.}

\item{full_data}{When False, only return the id and name fields (their specific names differ in specific nodes) for a node.
This is useful if you want your queries to return results faster with smaller amount of data requested.}

\item{limit}{Max number of items to retrieve.}

\item{offset}{Number of items to skip. Use \code{limit} and \code{offset} in combination to do pagination.}

\item{mode}{If \code{mode = "table"}, returns a data frame
(a \href{https://tibble.tidyverse.org/}{\code{tibble}} as per
\href{https://style.tidyverse.org/}{\code{tidyverse}} convention).
If \code{mode = "raw"}, returns a raw response from EpiGraphDB API
with minimal parsing done by \href{https://httr.r-lib.org/}{\code{httr}}.}
}
\value{
Data from \code{GET /meta/nodes/{meta_node}/list}
}
\description{
\href{https://docs.epigraphdb.org/api/api-endpoints/#get-metanodesmeta_nodelist}{\code{GET /meta/nodes/{meta_node}/list}}
}
\examples{
# List the first 5 Gwas nodes, with only id and name fields
\dontrun{
meta_nodes_list_node(meta_node = "Gwas", full_data = FALSE, limit = 5)
}

# List the 6th - 10th Disease nodes, with full properties
\dontrun{
meta_nodes_list_node(meta_node = "Disease", full_data = TRUE, limit = 5, offset = 0)
}
}
