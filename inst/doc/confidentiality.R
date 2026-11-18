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

## ----design-alias---------------------------------------------------------------------------------
df <- data.frame(
  site = factor(rep(c("Roseworthy", "Minnipa", "Turretfield"), each = 6L)),
  rep = rep(seq_len(3L), 6L),
  yield = rnorm(18L)
)
plan <- propose_roles(df, detect = FALSE)
plan <- set_role(plan, "site", role = "design", action = "alias")
s <- synthetic(mask(df, plan, mode = "collaborate", seed = 1L))
site_names_gone <- !any(as.character(s$site) %in% as.character(df$site))
site_names_gone
knitr::kable(
  as.data.frame(table(s$site), stringsAsFactors = FALSE),
  col.names = c("aliased site", "plots"),
  caption = "The three-site structure survives even though the labels do not."
)

## ----numeric-treatment----------------------------------------------------------------------------
rates <- data.frame(
  block = factor(rep(seq_len(4L), each = 5L)),
  nitro = rep(c(0, 0.2, 0.4, 0.6, 0.8), 4L)
)
rates$yield <- round(80 + 40 * rates$nitro + rnorm(20L, sd = 5), 1)
plan_rates <- propose_roles(rates, mode = "collaborate")
m_rates <- mask(rates, plan_rates, mode = "collaborate", seed = 1L)
knitr::kable(
  head(synthetic(m_rates), 5L),
  digits = 1,
  caption = "The first block of the clone: each N rate is an opaque label."
)
identical(unmask(synthetic(m_rates), recipe(m_rates))$nitro, rates$nitro)

## ----conditional-fixture--------------------------------------------------------------------------
set.seed(42)
n <- 600L
arm <- factor(rep(c("ctrl", "treat"), each = n / 2L))
# A real +5-unit effect of the treated arm on yield.
yield <- 10 + 5 * (arm == "treat") + rnorm(n, sd = 2)
trial <- data.frame(genotype = arm, yield = yield)

plan <- propose_roles(trial, detect = FALSE)
plan <- set_role(plan, "genotype", role = "treatment", action = "scramble")
plan <- set_role(plan, "yield", role = "outcome")

true_effect <- coef(lm(yield ~ genotype, trial))[["genotypetreat"]]
true_effect

## ----conditional-contrast-------------------------------------------------------------------------
marg <- mask(trial, plan, mode = "local", seed = 1L, conditional = FALSE)
cond <- mask(trial, plan, mode = "local", seed = 1L, conditional = TRUE)

effect_of <- function(m) {
  coef(lm(yield ~ genotype, synthetic(m)))[["genotypetreat"]]
}
effect_marg <- effect_of(marg)
effect_cond <- effect_of(cond)

knitr::kable(
  data.frame(
    clone = c("marginal (default)", "conditional"),
    estimated_effect = c(effect_marg, effect_cond),
    true_effect = true_effect
  ),
  digits = 3,
  caption = paste0(
    "The treated-arm effect, re-estimated on each clone against the ",
    "true effect fitted on the original trial."
  )
)

## ----conditional-marginal-held--------------------------------------------------------------------
knitr::kable(
  data.frame(
    source = c("original", "marginal clone", "conditional clone"),
    mean_yield = c(
      mean(trial$yield), mean(synthetic(marg)$yield),
      mean(synthetic(cond)$yield)
    ),
    sd_yield = c(
      sd(trial$yield), sd(synthetic(marg)$yield),
      sd(synthetic(cond)$yield)
    )
  ),
  digits = 3,
  caption = paste0(
    "Both clones preserve the pooled mean and SD of yield; only the ",
    "conditional clone preserves the arm-to-arm difference."
  )
)

## ----fig-conditional, eval = requireNamespace("ggplot2", quietly = TRUE), fig.cap = "Estimated treatment effect (yield units) on the marginal versus conditional clone, against the true effect fitted on the original trial (dashed line). The marginal clone collapses the effect toward zero; the conditional clone recovers it."----
effect_df <- data.frame(
  clone = factor(
    c("marginal (default)", "conditional"),
    levels = c("marginal (default)", "conditional")
  ),
  estimated_effect = c(effect_marg, effect_cond)
)
ggplot2::ggplot(
  effect_df, ggplot2::aes(x = clone, y = estimated_effect, fill = clone)
) +
  ggplot2::geom_col(width = 0.6) +
  ggplot2::geom_hline(yintercept = true_effect, linetype = "dashed") +
  ggplot2::scale_fill_viridis_d(guide = "none") +
  ggplot2::labs(
    x = "clone", y = "estimated treatment effect (yield units)",
    title = "Marginal versus conditional cloning of a known treatment effect"
  ) +
  ggplot2::theme_minimal()

## ----fig-conditional-skip, eval = !requireNamespace("ggplot2", quietly = TRUE), echo = FALSE, results = "asis"----
# cat("ggplot2 is not installed on this build, so the effect-comparison",
#     "figure is skipped; the table above carries the same numbers.\n")

## ----conditional-recipe---------------------------------------------------------------------------
recipe(cond)@conditional
recipe(cond)@conditioning_cols
recipe(cond)@conditioning_used
recipe(cond)@fallback_frac

## ----nonmonotone----------------------------------------------------------------------------------
set.seed(7)
n <- 2000L
x_lin <- rnorm(n)
y_lin <- 2 * x_lin + rnorm(n, sd = 0.5)     # monotone (linear)
x_u <- runif(n, -3, 3)
y_u <- x_u^2 + rnorm(n, sd = 0.5)           # non-monotone (U-shaped)
d <- data.frame(x_lin, y_lin, x_u, y_u)

s <- synthetic(mask(d, propose_roles(d, detect = FALSE), seed = 1L))

# R-squared of a quadratic fit captures association of any curvature.
fit_r2 <- function(x, y) summary(lm(y ~ poly(x, 2)))$r.squared
knitr::kable(
  data.frame(
    pair = c("monotone (y = 2x)", "non-monotone (y = x^2)"),
    original = c(fit_r2(x_lin, y_lin), fit_r2(x_u, y_u)),
    clone = c(fit_r2(s$x_lin, s$y_lin), fit_r2(s$x_u, s$y_u))
  ),
  digits = 3,
  caption = paste0(
    "R-squared of a quadratic fit, original versus clone, for a ",
    "monotone and a non-monotone pair."
  )
)

## ----grain----------------------------------------------------------------------------------------
set.seed(3)
trials <- data.frame(
  trial = factor(rep(sprintf("T%d", seq_len(8L)), each = 12L)),
  rep = factor(rep(rep(seq_len(3L), each = 4L), 8L)),
  gen = factor(rep(sprintf("G%d", seq_len(4L)), 24L))
)
trial_rain <- round(runif(8L, 250, 550))
trials$rainfall <- trial_rain[as.integer(trials$trial)]
trials$yield <- round(2 + trials$rainfall / 200 + rnorm(96L, sd = 0.4), 2)

plan_trials <- propose_roles(trials, mode = "collaborate")
by_unit <- mask(trials, plan_trials, mode = "collaborate", seed = 1L)
by_row <- mask(trials, plan_trials, mode = "collaborate", seed = 1L,
               grain = FALSE)
values_per_trial <- function(m) {
  s <- synthetic(m)
  as.vector(tapply(s$rainfall, s$trial, function(x) length(unique(x))))
}
knitr::kable(
  data.frame(
    clone = c("grain = TRUE (default)", "grain = FALSE"),
    distinct_rainfall_per_trial = c(
      paste(values_per_trial(by_unit), collapse = " "),
      paste(values_per_trial(by_row), collapse = " ")
    )
  ),
  caption = "Distinct rainfall values in each of the eight trials of the clone."
)

## ----events---------------------------------------------------------------------------------------
set.seed(5)
n_pat <- 300L
arm <- rep(c("A", "B"), length.out = n_pat)
t_rec <- rexp(n_pat, ifelse(arm == "A", 1 / 900, 1 / 1400))
t_death <- rexp(n_pat, 1 / 2500)
censor <- runif(n_pat, 30, 3000)
status_death <- as.integer(t_death <= censor)
time_death <- round(pmin(t_death, censor))
status_rec <- as.integer(t_rec < pmin(t_death, censor))
time_rec <- ifelse(status_rec == 1L, pmin(round(t_rec), time_death),
                   time_death)
patients <- data.frame(arm, time_rec, status_rec, time_death, status_death)

plan_pat <- propose_roles(patients, mode = "collaborate", detect = FALSE)
plan_pat <- set_role(plan_pat, "arm", role = "treatment")
for (cn in c("time_rec", "status_rec", "time_death", "status_death")) {
  plan_pat <- set_role(plan_pat, cn, role = "outcome", action = "scramble")
}
m_pat <- mask(
  patients, plan_pat, mode = "collaborate", seed = 1L,
  events = list(time = c("time_rec", "time_death"),
                status = c("status_rec", "status_death"))
)
s_pat <- synthetic(m_pat)
all(s_pat$time_rec <= s_pat$time_death)
pattern <- function(d) paste(d$status_rec, d$status_death, sep = "/")
knitr::kable(
  data.frame(
    pattern = names(table(pattern(patients))),
    original = as.vector(table(pattern(patients))),
    clone = as.vector(table(pattern(s_pat)))
  ),
  caption = "Rows in each recurrence/death status pattern."
)

## ----leaky-fixture--------------------------------------------------------------------------------
df <- data.frame(
  rep = rep(seq_len(3L), each = 20L),
  variety = factor(rep(paste0("V", seq_len(6L)), 10L)),
  contact_email = factor(rep(c("a@x", "b@y"), 30L)),
  rare_treatment = factor(c(
    "only_one",
    sample(c("alpha", "beta", "gamma"), 59L, replace = TRUE)
  )),
  yield = rnorm(60L, 5, 1),
  stringsAsFactors = FALSE
)

plan <- propose_roles(df, mode = "collaborate")
knitr::kable(
  plan[, c("col", "role", "action", "pii_suspected")],
  caption = paste0(
    "The proposed plan: `contact_email` is auto-flagged ",
    "`pii_suspected` and dropped."
  )
)

## ----leaky-mask-----------------------------------------------------------------------------------
plan <- set_role(plan, "yield", role = "outcome")
plan <- set_role(plan, "contact_email", role = "covariate", action = "keep")
plan <- set_role(plan, "rare_treatment", role = "covariate")
m <- mask(df, plan, mode = "collaborate", seed = 1L)
knitr::kable(
  audit_mask(m)[, c("col", "leakage_class", "notes")],
  caption = paste0(
    "The mask-time audit: `contact_email` and `rare_treatment` both ",
    "flag HIGH."
  )
)

## ----leaky-write, error = TRUE--------------------------------------------------------------------
try({
out_dir <- tempfile()
masque(df, roles = plan, out = out_dir, mode = "collaborate", seed = 1L,
       quiet = TRUE)
dir.exists(out_dir)
})

## ----kept-outcome---------------------------------------------------------------------------------
trial_kept <- data.frame(
  block = factor(rep(seq_len(4L), each = 6L)),
  variety = factor(rep(sprintf("V%d", seq_len(6L)), 4L)),
  height = round(rnorm(24L, 90, 8)),
  yield = round(rnorm(24L, 5, 1), 2)
)
plan_kept <- propose_roles(trial_kept, mode = "collaborate")
plan_kept <- set_role(plan_kept, "yield", role = "outcome", action = "keep")
m_kept <- mask(trial_kept, plan_kept, mode = "collaborate", seed = 1L)
knitr::kable(
  audit_mask(m_kept, print = FALSE)[col == "yield", c("col", "leakage_class")],
  caption = "A kept outcome in a collaborate clone."
)

## ----kept-outcome-allowed-------------------------------------------------------------------------
m_allowed <- mask(trial_kept, plan_kept, mode = "collaborate", seed = 1L,
                  allow_kept_outcome = TRUE)
knitr::kable(
  audit_mask(m_allowed, print = FALSE)[
    col == "yield", c("col", "leakage_class")
  ],
  caption = "The same outcome, allowed by the custodian."
)
recipe(m_allowed)@allow_kept_outcome

## ----small-table, error = TRUE--------------------------------------------------------------------
try({
tiny <- trial_kept[1:4, c("block", "yield")]
mask(tiny, propose_roles(tiny, mode = "collaborate"), mode = "collaborate",
     seed = 1L)
})

## ----set------------------------------------------------------------------------------------------
set_dir <- system.file("extdata", "met_set", package = "masque")
tabs <- read_set(set_dir)
tabs$agronomy[, year := factor(year)]
ms <- mask_set(tabs, mode = "collaborate", seed = 1L, quiet = TRUE)
ag <- synthetic(ms)$agronomy
qa <- synthetic(ms)$quality
join_survives <- setequal(unique(ag$gen), unique(qa$gen))
join_survives
head(unique(ag$gen), 3L)

## ----geo-coords, eval = requireNamespace("maps", quietly = TRUE)----------------------------------
sites <- data.frame(
  site = c("A", "B", "C"),
  lat  = c(-34.9, -35.2, -33.6),
  lon  = c(138.6, 142.0, 148.2)
)
knitr::kable(
  jitter_coordinates(sites, "lat", "lon", min_km = 5, max_km = 20, seed = 1L),
  digits = 4,
  caption = paste0(
    "Three sites, each displaced 5-20 km in a random direction and ",
    "re-drawn until on land."
  )
)

## ----geo-coords-skip, eval = !requireNamespace("maps", quietly = TRUE), echo = FALSE, results = "asis"----
# cat("maps is not installed on this build, so the coordinate-jitter demo",
#     "is skipped.\n")

## ----geo-by, eval = requireNamespace("maps", quietly = TRUE)--------------------------------------
plots <- data.frame(
  trial = rep(c("T1", "T2"), each = 4L),
  lat   = rep(c(-34.9, -35.2), each = 4L),
  lon   = rep(c(138.6, 142.0), each = 4L)
)
masked <- jitter_coordinates(plots, "lat", "lon", by = "trial", seed = 1L)
one_coord_per_trial <- nrow(unique(masked[c("lat", "lon")]))
one_coord_per_trial

## ----geo-by-skip, eval = !requireNamespace("maps", quietly = TRUE), echo = FALSE, results = "asis"----
# cat("maps is not installed on this build, so the site-grouping demo is",
#     "skipped.\n")

## ----geo-mask, eval = FALSE-----------------------------------------------------------------------
# plan <- set_role(plan, c("GPS_S", "GPS_E"), action = "keep")
# m <- mask(df, plan, mode = "collaborate",
#           coords = list(list(lat = "GPS_S", lon = "GPS_E", by = "TRIAL")),
#           seed = 1L)

## ----session-info---------------------------------------------------------------------------------
sessionInfo()

## ----reset-options, include = FALSE-------------------------------------------
options(old_options)

