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
has_agridat <- requireNamespace("agridat", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

## ----library--------------------------------------------------------------------------------------
library(masque)

## ----specs, eval = has_agridat--------------------------------------------------------------------
design_specs <- function() {
  augmented <- agridat::kling.augmented[, c("block", "gen", "tsw", "row", "col")]
  n_gen <- table(augmented$gen)
  augmented$check <- factor(ifelse(
    n_gen[as.character(augmented$gen)] > 1L, as.character(augmented$gen), "new"
  ))
  rcbd <- agridat::besag.elbatan
  rcbd$block <- factor(rcbd$col)
  list(
    CRD = list(
      data = agridat::cochran.crd, outcome = "inf", treatment = "trt",
      design = c("row", "col"),
      formula = inf ~ trt, terms = c(trt = "treatment")
    ),
    RCBD = list(
      data = rcbd, outcome = "yield", treatment = "gen",
      design = c("block", "row", "col"),
      formula = yield ~ block + gen,
      terms = c(block = "block", gen = "treatment")
    ),
    `Latin square` = list(
      data = agridat::fisher.latin, outcome = "yield", treatment = "trt",
      design = c("row", "col"),
      formula = yield ~ factor(row) + factor(col) + trt,
      terms = c(`factor(row)` = "block", `factor(col)` = "block",
        trt = "treatment")
    ),
    `split-plot` = list(
      data = agridat::yates.oats[, c("block", "row", "col", "gen", "nitro", "yield")],
      outcome = "yield", treatment = c("gen", "nitro"),
      design = c("block", "row", "col"),
      formula = yield ~ block + gen * nitro,
      terms = c(block = "block", gen = "treatment", nitro = "treatment")
    ),
    `alpha lattice` = list(
      data = agridat::john.alpha[, c("rep", "block", "gen", "yield", "row", "col")],
      outcome = "yield", treatment = "gen",
      design = c("rep", "block", "row", "col"),
      formula = yield ~ rep + rep:block + gen,
      terms = c(rep = "block", `rep:block` = "block", gen = "treatment")
    ),
    augmented = list(
      data = augmented, outcome = "tsw", treatment = c("gen", "check"),
      design = c("block", "row", "col"),
      formula = tsw ~ block + check,
      terms = c(block = "block", check = "treatment")
    ),
    MET = list(
      data = agridat::besag.met, outcome = "yield", treatment = "gen",
      design = c("county", "rep", "block", "row", "col"),
      formula = yield ~ county + county:rep + gen,
      terms = c(county = "environment", `county:rep` = "block",
        gen = "treatment")
    ),
    `repeated measures` = list(
      data = agridat::kenward.cattle, outcome = "weight", treatment = "trt",
      design = c("animal", "day"),
      formula = weight ~ trt + animal + factor(day),
      terms = c(trt = "treatment", animal = "block", `factor(day)` = "block")
    )
  )
}

## ----runner, eval = has_agridat-------------------------------------------------------------------
design_roles <- function(spec) {
  r <- propose_roles(spec$data, mode = "local", detect = FALSE)
  r <- set_role(r, spec$outcome, role = "outcome")
  r <- set_role(r, spec$treatment, role = "treatment", action = "keep")
  set_role(r, spec$design, role = "design", action = "keep")
}

level_means <- function(d, col, outcome) {
  tab <- table(d[[col]])
  m <- tapply(d[[outcome]], d[[col]], mean, na.rm = TRUE)
  m[names(tab)[tab >= 2L]]
}

design_run <- function(spec, seeds, ladder = "hierarchy") {
  d <- spec$data
  r <- design_roles(spec)
  f_orig <- anova(lm(spec$formula, data = d))
  terms <- names(spec$terms)
  alloc <- c(spec$treatment, spec$design)
  per_seed <- lapply(seeds, function(s) {
    m <- withCallingHandlers(
      mask(d, r, mode = "local", seed = s, conditional = TRUE, ladder = ladder),
      masque_conditional_degraded = function(w) invokeRestart("muffleWarning")
    )
    sy <- as.data.frame(synthetic(m))
    f_clone <- anova(lm(spec$formula, data = sy))
    trt <- t(vapply(spec$treatment, function(tc) {
      mo <- level_means(d, tc, spec$outcome)
      ms <- level_means(sy, tc, spec$outcome)[names(mo)]
      c(
        corr = if (length(mo) > 2L) cor(mo, ms) else NA_real_,
        ratio = sd(ms) / sd(mo),
        order = as.numeric(identical(names(which.max(mo)), names(which.max(ms))))
      )
    }, numeric(3)))
    list(
      alloc = all(vapply(alloc, function(cn) identical(sy[[cn]], d[[cn]]), logical(1))),
      na = identical(is.na(sy[[spec$outcome]]), is.na(d[[spec$outcome]])),
      F = setNames(f_clone[terms, "F value"], terms),
      trt = trt,
      used = paste(recipe(m)@conditioning_used, collapse = " x "),
      shifted = paste(recipe(m)@conditioning_shifted, collapse = ", ")
    )
  })
  F_clone <- do.call(rbind, lapply(per_seed, `[[`, "F"))
  trt <- simplify2array(lapply(per_seed, `[[`, "trt"))
  list(
    n = nrow(d),
    F_orig = setNames(f_orig[terms, "F value"], terms),
    F_clone = apply(F_clone, 2, median),
    alloc = all(vapply(per_seed, `[[`, logical(1), "alloc")),
    na = all(vapply(per_seed, `[[`, logical(1), "na")),
    trt_corr = apply(trt[, "corr", , drop = FALSE], 1, median),
    trt_ratio = apply(trt[, "ratio", , drop = FALSE], 1, median),
    trt_order = apply(trt[, "order", , drop = FALSE], 1, mean),
    used = per_seed[[1]]$used,
    shifted = per_seed[[1]]$shifted
  )
}

design_verdict <- function(spec, res, f_floor = 0.5, corr_floor = 0.7) {
  kind <- spec$terms
  block_terms <- names(kind)[kind != "treatment"]
  block_ok <- vapply(block_terms, function(tm) {
    res$F_orig[[tm]] < 2 || res$F_clone[[tm]] >= f_floor * res$F_orig[[tm]]
  }, logical(1))
  trt_detected <- any(res$F_orig[names(kind)[kind == "treatment"]] >= 2)
  trt_ok <- vapply(seq_along(spec$treatment), function(i) {
    corr <- res$trt_corr[[i]]
    signal <- if (is.na(corr)) res$trt_order[[i]] > 0.5 else corr >= corr_floor
    (signal || !trt_detected) && res$trt_ratio[[i]] >= 0.5 && res$trt_ratio[[i]] <= 2
  }, logical(1))
  c(allocation = res$alloc, na_mask = res$na,
    blocking = all(block_ok), treatment = all(trt_ok))
}

## ----run, eval = has_agridat----------------------------------------------------------------------
specs <- design_specs()
seeds <- 1:20
runs <- lapply(specs, design_run, seeds = seeds)
runs_levels <- lapply(specs, design_run, seeds = seeds, ladder = "levels")

## ----table, eval = has_agridat--------------------------------------------------------------------
verdicts <- t(mapply(design_verdict, specs, runs))
summary_tbl <- data.frame(
  design = names(specs),
  rows = vapply(runs, `[[`, numeric(1), "n"),
  stratum = vapply(runs, `[[`, character(1), "used"),
  carried = vapply(runs, `[[`, character(1), "shifted"),
  allocation = verdicts[, "allocation"],
  na_mask = verdicts[, "na_mask"],
  blocking = verdicts[, "blocking"],
  treatment = verdicts[, "treatment"],
  verdict = ifelse(rowSums(!verdicts) == 0, "pass", "fail"),
  row.names = NULL
)
knitr::kable(summary_tbl)

## ----terms, eval = has_agridat--------------------------------------------------------------------
term_tbl <- do.call(rbind, lapply(names(specs), function(nm) {
  data.frame(
    design = nm,
    term = names(specs[[nm]]$terms),
    kind = unname(specs[[nm]]$terms),
    F_original = round(runs[[nm]]$F_orig, 2),
    F_levels = round(runs_levels[[nm]]$F_clone, 2),
    F_hierarchy = round(runs[[nm]]$F_clone, 2),
    row.names = NULL
  )
}))
knitr::kable(term_tbl)

## ----treatment, eval = has_agridat----------------------------------------------------------------
trt_tbl <- do.call(rbind, lapply(names(specs), function(nm) {
  data.frame(
    design = nm,
    treatment = specs[[nm]]$treatment,
    corr_levels = round(runs_levels[[nm]]$trt_corr, 2),
    corr_hierarchy = round(runs[[nm]]$trt_corr, 2),
    spread_ratio = round(runs[[nm]]$trt_ratio, 2),
    row.names = NULL
  )
}))
knitr::kable(trt_tbl)

## ----figure, eval = has_agridat && has_ggplot2, fig.height = 5.5, fig.cap = "F statistic of each model term on the clone as a share of the original, median over twenty seeds, by ladder. Points on the dashed line are terms the clone reproduces exactly; the dotted line is the pass floor of one half. The levels ladder loses most of every blocking and environment term; the hierarchy ladder keeps them.", fig.alt = "Dot plot of clone-to-original F ratios for every design term, two ladders side by side."----
plot_df <- rbind(
  data.frame(term_tbl[, c("design", "term", "kind")], ladder = "levels",
    ratio = term_tbl$F_levels / term_tbl$F_original),
  data.frame(term_tbl[, c("design", "term", "kind")], ladder = "hierarchy",
    ratio = term_tbl$F_hierarchy / term_tbl$F_original)
)
plot_df <- plot_df[term_tbl$F_original >= 2, ]
plot_df$label <- paste(plot_df$design, plot_df$term, sep = ": ")
ggplot2::ggplot(plot_df, ggplot2::aes(x = ratio, y = label, colour = ladder)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey50") +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = c(levels = "#E69F00", hierarchy = "#0072B2")) +
  ggplot2::facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
  ggplot2::labs(x = "F on the clone / F on the original (log scale)", y = NULL) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text.y = ggplot2::element_text(angle = 0)
  )

## ----no-agridat, eval = !has_agridat, echo = FALSE, results = "asis"------------------------------
# cat("The agridat package is not installed on this build, so the designs",
#     "are not masked here; the tables and the figure are omitted.\n")

## ----read-values, eval = has_agridat, include = FALSE---------------------------------------------
n_pass <- sum(summary_tbl$verdict == "pass")
met <- runs$MET
met_levels <- runs_levels$MET
oats <- runs$`split-plot`
oats_levels <- runs_levels$`split-plot`

## ----session-info---------------------------------------------------------------------------------
sessionInfo()

## ----reset-options, include = FALSE-------------------------------------------
options(old_options)

