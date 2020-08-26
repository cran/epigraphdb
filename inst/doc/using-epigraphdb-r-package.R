## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("epigraphdb")

## -----------------------------------------------------------------------------
df <- mr(
  exposure_trait = "Body mass index",
  outcome_trait = "Coronary heart disease",
  mode = "table"
)
df

## -----------------------------------------------------------------------------
df <- query_epigraphdb(
  route = "/mr",
  params = list(
    exposure_trait = "Body mass index",
    outcome_trait = "Coronary heart disease"
  ),
  mode = "table"
)

df

## -----------------------------------------------------------------------------
df <- mr(
  exposure_trait = "Body mass index",
  outcome_trait = "Coronary heart disease"
)
df

## -----------------------------------------------------------------------------
response <- mr(
  exposure_trait = "Body mass index",
  outcome_trait = "Coronary heart disease",
  mode = "raw"
)
response %>% str(2)

