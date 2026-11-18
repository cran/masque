## ----setup, include = FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150
)
old_options <- options(width = 100)
set.seed(1)

## ----library--------------------------------------------------------------------------------------
library(masque)

## ----one-call-------------------------------------------------------------------------------------
f <- system.file("extdata", "john_alpha.csv", package = "masque")
df <- read.csv(f, stringsAsFactors = TRUE)
knitr::kable(
  head(df),
  caption = "The first six plots of the alpha-design trial."
)

m <- masque(df, mode = "collaborate", seed = 1L, ask = FALSE)

## ----one-call-extract-----------------------------------------------------------------------------
synth <- synthetic(m)
knitr::kable(
  head(synth),
  caption = paste0(
    "The first six plots of the synthetic clone: `gen` is aliased, ",
    "`yield` is re-simulated."
  )
)

## ----plan-----------------------------------------------------------------------------------------
plan <- propose_roles(df, mode = "collaborate")
knitr::kable(
  plan[, c("col", "role", "action", "kind")],
  caption = paste0(
    "The proposed masking plan for the alpha-design trial ",
    "(collaborate mode)."
  )
)

## ----set-role-------------------------------------------------------------------------------------
plan_id <- set_role(plan, "plot", role = "id")
knitr::kable(
  plan_id[, c("col", "role", "action")],
  caption = "The plan after re-roling `plot` as an identifier."
)

## ----edit-vectorised------------------------------------------------------------------------------
r2 <- set_role(plan, c("row", "col"), action = "drop")

## ----edit-direct----------------------------------------------------------------------------------
r2[col == "rep", `:=`(role = "covariate", action = NA_character_)]
knitr::kable(
  r2[, c("col", "role", "action")],
  caption = paste0(
    "A hand-edited plan: `row`/`col` dropped, `rep` re-roled to a ",
    "covariate."
  )
)

## ----role-options---------------------------------------------------------------------------------
knitr::kable(
  role_options(kind = "factor"),
  caption = paste0(
    "Every (role, action) pair the validator accepts for a factor ",
    "column."
  )
)

## ----pii-by-hand, eval = FALSE--------------------------------------------------------------------
# plan[col == "comments", pii_suspected := TRUE]

## ----mask-----------------------------------------------------------------------------------------
m <- mask(df, plan, mode = "collaborate", seed = 1L)
synth <- synthetic(m)
byte_identical_plot <- identical(synth$plot, df$plot)
gen_levels_changed <- !setequal(levels(synth$gen), levels(df$gen))
byte_identical_plot
gen_levels_changed
head(levels(synth$gen))

## ----fail-closed, error = TRUE--------------------------------------------------------------------
try({
bad_plan <- set_role(plan, "plot", role = "design", action = "scramble")
mask(df, bad_plan, mode = "collaborate", seed = 1L)
})

## ----set------------------------------------------------------------------------------------------
set_dir <- system.file("extdata", "met_set", package = "masque")
ms <- masque(set_dir, mode = "collaborate", seed = 1L, ask = FALSE)
ms

## ----set-year-------------------------------------------------------------------------------------
tabs <- read_set(set_dir)
tabs$agronomy[, year := factor(year)]
ms <- masque(tabs, mode = "collaborate", seed = 1L, ask = FALSE, quiet = TRUE)
knitr::kable(
  head(synthetic(ms)$agronomy[, c("env", "loc", "year", "gen")], 3L),
  caption = "The first rows of the clone, with `year` now aliased."
)

## ----data-card-roles------------------------------------------------------------------------------
knitr::kable(
  roles(m),
  caption = "The role table of the clone, as the data card reads it."
)

## ----data-card------------------------------------------------------------------------------------
card_file <- tempfile(fileext = ".yaml")
card <- data_card(
  m,
  file = card_file,
  name = "alpha_trial",
  spatial = c("row", "col"),
  unit = "plot"
)
writeLines(head(readLines(card_file), 30))

## ----met-scope------------------------------------------------------------------------------------
met <- expand.grid(
  env = factor(c("E1", "E2")),
  rep = factor(seq_len(2L)),
  gen = factor(c("G1", "G2", "G3"))
)
met$yield <- seq_len(nrow(met)) + rep(c(0, 2), each = 6L)

ds <- detect_design(met)
ds@scope_label
ds@env_cols
ds@connectivity$status
ds@within_design_label

## ----met-explicit---------------------------------------------------------------------------------
ds_explicit <- detect_design(met, env = "env")

## ----met-plot, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Genotypes observed per environment (E1, E2) for the toy multi-environment trial; the subtitle carries connectivity status, component count, and the recovered within-environment design."----
plot(ds_explicit, df = met, engine = "ggplot2")

## ----met-plot-skip, eval = !requireNamespace("ggplot2", quietly = TRUE), echo = FALSE, results = "asis"----
# cat("ggplot2 is not installed on this build, so the coverage figure is",
#     "skipped; `plot(ds_explicit, df = met)` still draws the base-graphics",
#     "equivalent.\n")

## ----met-roles------------------------------------------------------------------------------------
met_roles <- propose_roles(met, mode = "collaborate")
knitr::kable(
  met_roles[met_roles$col == "env", c("col", "role", "action")],
  caption = "How the masking plan treats the `env` column."
)

## ----audit-tally, include = FALSE-----------------------------------------------------------------
audit <- audit_mask(m, print = FALSE)
grades <- c(audit$leakage_class, attr(audit, "joint")$leakage_class)
n_high <- sum(grades == "high")
n_medium <- sum(grades == "medium")

## ----session-info---------------------------------------------------------------------------------
sessionInfo()

## ----reset-options, include = FALSE-------------------------------------------
options(old_options)

