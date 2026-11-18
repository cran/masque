## ----setup, include = FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150
)
old_options <- options(width = 100)
library(masque)
set.seed(1)

f <- system.file("extdata", "john_alpha.csv", package = "masque")
df <- read.csv(f, stringsAsFactors = TRUE)
plan <- propose_roles(df, mode = "collaborate")

## ----accessors------------------------------------------------------------------------------------
m <- mask(df, plan, mode = "collaborate", seed = 1L)
rec <- recipe(m)
class(rec)

## ----coords-account, eval = requireNamespace("maps", quietly = TRUE)------------------------------
df_geo <- df
df_geo$lat <- -34.9 + stats::runif(nrow(df_geo), -0.05, 0.05)
df_geo$lon <- 138.6 + stats::runif(nrow(df_geo), -0.05, 0.05)
plan_geo <- propose_roles(df_geo, mode = "collaborate")
plan_geo <- set_role(plan_geo, c("lat", "lon"), action = "keep")

m_geo <- mask(df_geo, plan_geo, mode = "collaborate", seed = 1L,
              coords = list(list(lat = "lat", lon = "lon")))
n_sites_recorded <- recipe(m_geo)@coords[[1]]$n_sites
recipe(m_geo)@coords[[1]][c("method", "min_km", "max_km", "n_sites")]
recipe(m_geo)@allow_unmasked_coords

## ----coords-account-skip, eval = !requireNamespace("maps", quietly = TRUE), echo = FALSE, results = "asis"----
# cat("maps is not installed on this build, so the coordinate-account demo",
#     "is skipped.\n")
# n_sites_recorded <- NA_integer_

## ----print-redacted-------------------------------------------------------------------------------
rec

## ----reveal, results = "hide"---------------------------------------------------------------------
reveal_maps(rec)

## ----roundtrip-pipeline---------------------------------------------------------------------------
# Analyst trains on the synthetic clone.
fit <- lm(yield ~ gen, data = synthetic(m))

# Custodian re-targets the same pipeline onto the real data.
orig_in_synth_space <- apply_recipe(df, rec)
preds <- predict(fit, newdata = orig_in_synth_space)
preds_match_nrow <- length(preds) == nrow(df)
preds_match_nrow

## ----roundtrip-labels-----------------------------------------------------------------------------
fwd <- apply_recipe(df, rec)
back <- unmask(fwd, rec)
labels_recovered <- identical(as.character(back$gen), as.character(df$gen))
labels_recovered

## ----fail-closed, error = TRUE--------------------------------------------------------------------
try({
drifted <- df
levels(drifted$gen) <- c(levels(df$gen), "BRAND_NEW")
drifted$gen[1] <- "BRAND_NEW"
apply_recipe(drifted, rec)
})

## ----save-----------------------------------------------------------------------------------------
tmp <- tempfile(fileext = ".rds")
save_recipe(rec, tmp)
rec2 <- read_recipe(tmp)
fingerprint_survives_roundtrip <- identical(rec@integrity_fp, rec2@integrity_fp)
fingerprint_survives_roundtrip

## ----set-bundle-----------------------------------------------------------------------------------
set_dir <- system.file("extdata", "met_set", package = "masque")
originals <- read_set(set_dir)
originals$agronomy[, year := factor(year)]
ms <- mask_set(originals, mode = "collaborate", seed = 1L, quiet = TRUE)
fwd_set <- apply_recipe(originals, recipe(ms))
names(fwd_set)

## ----fig-overlay, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Density of yield, original trial against the synthetic clone the analyst develops against. The two distributions overlap closely because the default numeric synthesis preserves each column's marginal distribution."----
synth_tbl <- synthetic(m)
overlay_df <- rbind(
  data.frame(value = df$yield, table = "original"),
  data.frame(value = synth_tbl$yield, table = "synthetic")
)
ggplot2::ggplot(overlay_df, ggplot2::aes(x = value, fill = table)) +
  ggplot2::geom_density(alpha = 0.5) +
  ggplot2::scale_fill_viridis_d(name = "table") +
  ggplot2::labs(
    x = "yield (t/ha)", y = "density",
    title = "Original versus synthetic marginal distribution of yield"
  ) +
  ggplot2::theme_minimal()

## ----fig-overlay-skip, eval = !requireNamespace("ggplot2", quietly = TRUE), echo = FALSE, results = "asis"----
# cat("ggplot2 is not installed on this build, so the marginal-overlay",
#     "figure is skipped.\n")

## ----session-info---------------------------------------------------------------------------------
sessionInfo()

## ----reset-options, include = FALSE-------------------------------------------
options(old_options)

