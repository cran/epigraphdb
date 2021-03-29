## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("dplyr")
library("epigraphdb")

## -----------------------------------------------------------------------------
# Let's test this on a few genes
genes <- c("TP53", "BRCA1", "TNF")

# Select the endpoint "POST /mappings/gene-to-protein"
endpoint <- "/mappings/gene-to-protein"
method <- "POST"

# Build a query
proteins <- query_epigraphdb(
  route = endpoint,
  params = list(gene_name_list = genes),
  mode = "table",
  method = method
)

# Check the output
print(proteins)

## -----------------------------------------------------------------------------
# Let's see what pathways the protein product of TP53 gene is involved in

# NOTE: Argument `uniprot_id_list` requires a list of UniProt IDs in
#       POST /mappings/gene-to-protein.
#       In this case for the R package, we need to wrap `proteins_uniprot_ids`
#       with an `I()` function (AsIs) to prevent auto-unpacking by `httr`
proteins_uniprot_ids <- c("P04637") %>% I()

endpoint <- "/protein/in-pathway"

pathway_df <- query_epigraphdb(
  route = endpoint,
  params = list(uniprot_id_list = proteins_uniprot_ids),
  mode = "table",
  method = "POST"
)

# Check out how many pathways this protein is found in
print(pathway_df)

# Get pathways names (Reactome IDs)
print(pathway_df$pathway_reactome_id[[1]])

## -----------------------------------------------------------------------------
# Let's see what exactly this pathway is
reactome_id <- "R-HSA-6804754"

endpoint <- "/meta/nodes/Pathway/search"

pathway_info <- query_epigraphdb(
  route = endpoint,
  params = list(id = reactome_id),
  mode = "table"
)

# Pathway description
print(pathway_info)

## -----------------------------------------------------------------------------
mappings_gene_to_protein(genes)
protein_in_pathway(proteins_uniprot_ids)

## -----------------------------------------------------------------------------
# Let's see what Body mass index GWAS are available in EpiGraphDB
trait <- "body mass index"

endpoint <- "/meta/nodes/Gwas/search"

results <- query_epigraphdb(
  route = endpoint,
  params = list(name = trait),
  mode = "table"
)

# show selected columns in the results
results %>%
  select(
    node.trait, node.id, node.sample_size,
    node.year, node.author
  )

## -----------------------------------------------------------------------------
# Look up all MR analyses where a trait was used as exposure
# and find all outcome traits with causal evidence from it.

trait1 <- "Body mass index"
# NB: here, trait name has to specific (not fuzzy):
# use the exact trait name wording as in GWAS `node.trait` (previous example)

endpoint <- "/mr"

mr_df <- query_epigraphdb(
  route = endpoint,
  params = list(
    exposure_trait = trait1,
    pval_threshold = 5e-08
  ),
  mode = "table"
)
print(mr_df)

## -----------------------------------------------------------------------------
# Show how many MR analyses were done for each Body mass index GWAS
mr_df %>% count(exposure.id)

## -----------------------------------------------------------------------------
# Look up all MR for a specified outcome trait
# and find all exposure traits with causal evidence on it.

trait2 <- "Waist circumference"

endpoint <- "/mr"

mr_df <- query_epigraphdb(
  route = endpoint,
  params = list(
    outcome_trait = trait2,
    pval_threshold = 5e-08
  ),
  mode = "table"
)
print(mr_df)

## -----------------------------------------------------------------------------
# Look up a specific pair of exposure+outcome
trait1 <- "Body mass index"
trait2 <- "Coronary heart disease"

endpoint <- "/mr"

mr_df <- query_epigraphdb(
  route = endpoint,
  params = list(
    exposure_trait = trait1,
    outcome_trait = trait2
  ),
  mode = "table"
)

print(mr_df)

## -----------------------------------------------------------------------------
mr(
  exposure_trait = trait1,
  outcome_trait = trait2
)

## -----------------------------------------------------------------------------
# Identity all publications where a gene is mentioned with relation to a disease

gene <- "IL23R"
trait <- "Inflammatory bowel disease"

endpoint <- "/gene/literature"

lit_df <- query_epigraphdb(
  route = endpoint,
  params = list(
    gene_name = gene,
    object_name = trait
  ),
  mode = "table"
)

# Review the found evidence in the literature
print(lit_df)

## -----------------------------------------------------------------------------
# Get a list of all PubMed IDs
lit_df %>%
  pull(pubmed_id) %>%
  unlist() %>%
  unique()

## -----------------------------------------------------------------------------
endpoint <- "/meta/nodes/list"
meta_node_fields <- query_epigraphdb(
  route = endpoint, params = NULL, mode = "raw"
)
meta_node_fields %>% unlist()

## -----------------------------------------------------------------------------
name <- "breast cancer"

endpoint <- "/meta/nodes/Gwas/search"
results <- query_epigraphdb(
  route = endpoint,
  params = list(name = name),
  mode = "table"
)
results %>%
  select(
    node.trait, node.id, node.sample_size,
    node.year, node.author
  )

endpoint <- "/meta/nodes/Disease/search"
results <- query_epigraphdb(
  route = endpoint,
  params = list(name = name),
  mode = "table"
)
results %>%
  select(node.label, node.id, node.definition)

endpoint <- "/meta/nodes/Efo/search"
results <- query_epigraphdb(
  route = endpoint,
  params = list(name = name),
  mode = "table"
)
results

## -----------------------------------------------------------------------------
# Running a GWAS query from Part 2
# Ask for "raw" format (as a list)
trait <- "body mass index"
endpoint <- "/meta/nodes/Gwas/search"
response <- query_epigraphdb(
  route = endpoint,
  params = list(name = trait),
  mode = "raw"
)

# display the Cypher query
response$metadata$query

## -----------------------------------------------------------------------------
query <- "
    MATCH (node: Gwas)
    WHERE node.trait =~ \"(?i).*body mass index.*\"
    RETURN node LIMIT 10;
"

## -----------------------------------------------------------------------------
# use POST /cypher

endpoint <- "/cypher"
method <- "POST"
params <- list(query = query)

results <- query_epigraphdb(
  route = endpoint,
  params = params,
  method = method,
  mode = "table"
)

# The result should be identical to the example in Part 2
results %>%
  select(
    node.trait, node.id, node.sample_size,
    node.year, node.author
  )

## -----------------------------------------------------------------------------
# Let's return Body mass index GWAS that were done by Locke AE

query <- "
    MATCH (node: Gwas)
    WHERE node.trait =~ \"(?i).*body mass index.*\"
    AND node.author = \"Locke AE\"
    RETURN node;
"

endpoint <- "/cypher"
method <- "POST"
params <- list(query = query)

results_subset <- query_epigraphdb(
  route = endpoint,
  params = params,
  method = method,
  mode = "table"
)

results_subset %>%
  select(
    node.trait, node.id, node.sample_size,
    node.year, node.author
  )

## -----------------------------------------------------------------------------
# Extract info only for GWAS 'ieu-a-2'
query <- "
    MATCH (node: Gwas)
    WHERE node.id = \"ieu-a-2\"
    RETURN node;
"
endpoint <- "/cypher"
params <- list(query = query)
results_subset <- query_epigraphdb(
  route = endpoint,
  params = params,
  method = "POST",
  mode = "table"
)

results_subset %>%
  select(
    node.trait, node.id, node.sample_size,
    node.year, node.author
  )

## -----------------------------------------------------------------------------
# Return MR results only for exposure trait 'ieu-a-2' (body mass index)

# first let's check the MR Cypher query that we run in "table" mode in Part 2
trait1 <- "Body mass index"
endpoint <- "/mr"
mr_df <- query_epigraphdb(
  route = endpoint,
  params = list(
    exposure_trait = trait1,
    pval_threshold = 5e-08
  ),
  mode = "raw"
)
mr_df$metadata$query

# modify the query to only return 'ieu-a-2' GWAS results
query <- "
  MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas)
  WHERE exposure.id = \"ieu-a-2\"
  AND mr.pval < 5e-08
  RETURN exposure {.id, .trait}, outcome {.id, .trait}, mr {.b, .se, .pval, .method, .selection, .moescore}
  ORDER BY mr.pval ;
"

endpoint <- "/cypher"
params <- list(query = query)
results_subset <- query_epigraphdb(
  route = endpoint,
  params = params,
  method = "POST",
  mode = "table"
)
results_subset

# check exposures in the results
results_subset %>% count(exposure.id)

