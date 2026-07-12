# =============================================================================
# Modeling methods library
# -----------------------------------------------------------------------------
# Every modeling method is defined here as a single documented function. A
# compute orchestrator (`compute_all_methods()`) applies them to a set of
# conditions; a separate writer (`write_all_methods()`) exports the CSVs.
# Forest plots follow the same split: `plot_all_methods()` builds the figures
# (selective LASSO at both lambdas + split), `write_all_plots()` saves them.
#
# Methods and their column prefixes in the combined table:
#   * lasso_opt_ / lasso_min_  selective-inference LASSO (fit_selective_lasso)
#   * boot_lasso_              bootstrap of the selective-inference LASSO
#   * split_                   multi-sample splitting with cross-fitting
#   * firth_                   Firth penalized logistic on the selected set
#   * step_ / fs_             forward stepwise (logistic refit / native linear)
#   * multi_                  p-value-screened multivariable logistic
#   * uni_                    univariable logistic
#
# Reading guide for the tricky parts:
#   - Selective inference: `fit_selective_lasso()` corrects p-values/CIs for the
#     fact that the LASSO chose the variables. Its comments explain the lambda
#     scaling, the standardization choice, and the stable-plateau rule.
#   - Split (`fit_split_lasso()`): the only method whose p-values are BOTH
#     selection-corrected AND on a scale a clinician reads directly, because it
#     selects on one half of the data and tests on the other (never the same
#     rows), then aggregates over many random splits (Meinshausen, Meier &
#     Buhlmann 2009). Cross-fitting uses each half for both roles.
# =============================================================================


# ---- Formatting helpers -----------------------------------------------------

#' Format an odds ratio (or CI bound) as a display string
#'
#' @description
#' One rule for every OR-scale number, so plots and tables always match:
#'   * `>= 1`: three significant figures, keeping trailing zeros (1.6 -> "1.60",
#'     129 -> "129", 12.5 -> "12.5").
#'   * `< 1`: two decimal places (0.905 -> "0.91", 0.808 -> "0.81"). A bound
#'     below 0.01 therefore reads "0.00", which is intended: at that magnitude
#'     the exact value carries no clinical meaning and the extra digits only
#'     lengthen the label.
#'   * `>= 1000`: scientific with two significant figures, to stay compact.
#' @param v Numeric vector.
#' @return Character vector; `NA` where `v` is `NA`.
#' @export
fmt_or_num <- function(v) {
  vapply(v, function(x) {
    if (is.na(x)) return(NA_character_)
    ax <- abs(x)
    if (ax >= 1000) return(formatC(x, format = "e", digits = 1))
    if (ax >= 1) {
      # Decimals needed for three significant figures at this magnitude
      return(formatC(x, format = "f", digits = 2 - floor(log10(ax))))
    }
    formatC(x, format = "f", digits = 2)
  }, character(1))
}


#' Numeric sibling of [fmt_or_num()] (same rounding, but returns a number)
#'
#' Use where a numeric column must stay numeric (e.g. before [ci_str()] rebuilds
#' the string); the display formatter [fmt_or_num()] applied later is then a
#' no-op, so there is no double rounding.
#' @param x Numeric vector.
#' @return Numeric vector, rounded to the OR display rule.
#' @export
round_or <- function(x) {
  ifelse(abs(x) < 1, round(x, 2), signif(x, 3))
}


#' Combine two bounds into one readable "lo - hi" string
#'
#' @param lo,hi Numeric vectors of lower and upper bounds.
#' @return Character vector; `NA` where either bound is `NA`.
#' @export
ci_str <- function(lo, hi) {
  # Return NA when a bound is missing, else the OR-formatted "lo - hi" string
  ifelse(
    is.na(lo) | is.na(hi),
    NA_character_,
    paste0(fmt_or_num(lo), " - ", fmt_or_num(hi))
  )
}


#' Prefix every non-`term` column of a block
#'
#' @param df A tibble with a `term` column.
#' @param pre Prefix string (e.g. `"lasso_opt"`).
#' @return `df` with columns renamed `<pre>_<col>`.
#' @export
prefix_cols <- function(df, pre) {
  # Rename all columns except `term` to `<pre>_<col>`
  dplyr::rename_with(df, ~ paste0(pre, "_", .x), -term)
}


#' Round all numeric columns to a readable number of significant digits
#'
#' @param df A tibble.
#' @param digits Significant digits (keeps small p-values legible).
#' @return `df` with numeric columns passed through `signif()`.
#' @export
round_readable <- function(df, digits = 3) {
  # Apply signif() to every numeric column (leaves strings/CIs untouched)
  dplyr::mutate(df, dplyr::across(dplyr::where(is.numeric), ~ signif(.x, digits)))
}


#' Turn raw model.matrix term names into readable labels
#'
#' @param old_terms Character vector of raw term names.
#' @param col_names Named lookup vector (name = label, value = variable).
#' @return Character vector of labels; unmatched terms fall back to the raw name.
#' @importFrom stringr str_extract str_replace
#' @importFrom purrr map_chr
#' @export
prettify_terms <- function(old_terms, col_names) {
  # Trailing capitalised part is the factor level (e.g. "TRUE", "Deep")
  factor_terms <- stringr::str_extract(old_terms, "[A-Z].*")
  # Base variable name is what remains after stripping that level
  base_terms <- stringr::str_replace(old_terms, "[A-Z].*$", "")
  # Look up the human label; fall back to the raw name when not found
  new_terms <-
    purrr::map_chr(base_terms, function(bt) {
      label <- names(col_names)[col_names == bt]
      if (length(label) == 0) bt else label[1]
    })
  # Re-attach the factor level when there was one
  ifelse(
    is.na(factor_terms),
    new_terms,
    paste(new_terms, factor_terms, sep = " - ")
  )
}


#' Normalize a glm/broom style OR table to the common schema
#'
#' @param tbl Tibble with columns `Predictors`, `Odds Ratios (OR)`, `SE`,
#'   `CI (low)`, `CI (high)`, `P-values` (the univariable/multivariable format).
#' @return Tibble of `term`, `OR`, `SE`, `ci`, `pval`.
#' @export
norm_orlike <- function(tbl) {
  # Map the stylised univariable/multivariable columns onto the shared schema
  tbl %>%
    dplyr::transmute(
      term = .data[["Predictors"]],
      OR = .data[["Odds Ratios (OR)"]],
      SE = .data[["SE"]],
      ci = ci_str(.data[["CI (low)"]], .data[["CI (high)"]]),
      pval = .data[["P-values"]]
    )
}


# ---- Internal fitters (small, shared building blocks) -----------------------

#' Select variables by LASSO at lambda.min (internal)
#'
#' @param x Centered design matrix.
#' @param y 0/1 outcome.
#' @param nfolds CV folds.
#' @param elastic_net_alpha glmnet mixing (1 = lasso).
#' @param standardize glmnet column standardization.
#' @return Character vector of selected variable names (no intercept).
#' @keywords internal
select_lasso_vars <- function(
    x,
    y,
    nfolds = 10,
    elastic_net_alpha = 1,
    standardize = FALSE
) {
  # Cross-validate to find lambda.min
  cv <-
    glmnet::cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      alpha = elastic_net_alpha,
      nfolds = nfolds,
      standardize = standardize
    )
  # Read the coefficients at lambda.min
  b <- as.matrix(coef(cv, s = "lambda.min"))
  # Return the non-zero, non-intercept variable names
  setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
}


#' Logistic estimate + Wald p-values on a fixed variable set (internal)
#'
#' @param x Design matrix.
#' @param y 0/1 outcome.
#' @param vars Variable names to include.
#' @return List of `terms`, `coef` (log-OR), `se`, `pval`; `NULL` on failure.
#' @keywords internal
glm_estimate <- function(x, y, vars) {
  # Nothing to estimate if no variables were passed
  if (length(vars) == 0) return(NULL)
  # Build a data frame with the outcome and the selected columns
  df <- data.frame(.y = y, x[, vars, drop = FALSE], check.names = FALSE)
  # Fit ordinary logistic regression (may warn under separation)
  fit <- tryCatch(
    suppressWarnings(glm(.y ~ ., data = df, family = binomial())),
    error = function(e) NULL
  )
  # Bail out gracefully if the fit failed
  if (is.null(fit)) return(NULL)
  # Pull the coefficient table and drop the intercept row
  co <- summary(fit)$coefficients
  keep <- setdiff(rownames(co), "(Intercept)")
  # Bail out if only the intercept survived
  if (length(keep) == 0) return(NULL)
  # Return log-ORs, SEs and p-values aligned to the kept terms
  list(
    terms = keep,
    coef = co[keep, "Estimate"],
    se = co[keep, "Std. Error"],
    pval = co[keep, "Pr(>|z|)"]
  )
}


#' Firth-penalized logistic estimate on a fixed variable set (internal)
#'
#' @description
#' Same contract as [glm_estimate()] but fits via `logistf::logistf()`, so
#' estimates stay finite under separation and small-sample p-values are better
#' calibrated. Used by the split method, whose half-samples are small and
#' separation-prone — exactly where ordinary MLE is unreliable.
#'
#' @param x Design matrix.
#' @param y 0/1 outcome.
#' @param vars Variable names to include.
#' @return List of `terms`, `coef` (log-OR), `se`, `pval` (penalized-LR);
#'   `NULL` on failure.
#' @keywords internal
firth_estimate <- function(x, y, vars) {
  # Nothing to estimate if no variables were passed
  if (length(vars) == 0) return(NULL)
  # Build a data frame with the outcome and the selected columns
  df <- data.frame(.y = y, x[, vars, drop = FALSE], check.names = FALSE)
  # Firth-penalized logistic (finite under separation); NULL on failure
  fit <- tryCatch(
    suppressWarnings(logistf::logistf(.y ~ ., data = df)),
    error = function(e) NULL
  )
  # Bail out gracefully if the fit failed
  if (is.null(fit)) return(NULL)
  # Non-intercept coefficient positions
  idx <- which(names(stats::coef(fit)) != "(Intercept)")
  # Bail out if only the intercept survived
  if (length(idx) == 0) return(NULL)
  # Return log-ORs, SEs and penalized-LR p-values aligned to the kept terms
  list(
    terms = names(stats::coef(fit))[idx],
    coef = stats::coef(fit)[idx],
    se = sqrt(diag(stats::vcov(fit)))[idx],
    pval = fit$prob[idx]
  )
}


# ---- Method 1: selective-inference LASSO ------------------------------------

#' Selective-inference LASSO for a binary outcome
#'
#' @description
#' Fit an L1-penalized logistic regression, choose a regularization-stable model
#' (or a user-specified one), and return post-selection (selective) p-values and
#' confidence intervals via [selectiveInference::fixedLassoInf()]. This is the
#' method that formally corrects the inference for the model-selection step.
#'
#' @details
#' Design decisions baked in from debugging on small, separation-prone data:
#'
#' * **Centered, unscaled `x`.** Scaling the columns (equivalently
#'   `glmnet(standardize = TRUE)`) pushed the observed statistics onto the
#'   boundary of the truncation region and produced infinite intervals.
#' * **`standardize = FALSE`**, so glmnet penalizes exactly the design handed to
#'   `fixedLassoInf` with a uniform penalty.
#' * **`lambda * n`.** `fixedLassoInf` uses the standard lasso objective while
#'   glmnet scales the loss by `1/n`; `beta` is read at glmnet's `lambda`, so the
#'   matching penalty is `lambda * n` (otherwise the KKT check fails).
#' * **Fixed folds** via an explicit `foldid` for reproducibility.
#' * **Stable-plateau selection.** A fractional-SE sweep generalises the 1-SE
#'   rule; with `k = NULL` the function takes the longest run of tolerances that
#'   select an identical variable set. `k = 0` reproduces `lambda.min`, `k = 1`
#'   reproduces `lambda.1se`. This is a pre-specified rule, not tuned to results.
#'
#' The reported `pval` is the two-sided selective p (`2 * min(p, 1 - p)`);
#' `pval_1sided` is the package-native one-sided p (one-sided because the LASSO
#' already fixed
#' the sign, so the test conditions on that direction).
#'
#' @param x Numeric design matrix (`n` x `p`), PRE-CENTERED, NOT column-scaled.
#' @param y Numeric 0/1 outcome of length `n`.
#' @param family glmnet family; validated for `"binomial"`.
#' @param elastic_net_alpha glmnet mixing (1 = lasso, 0 = ridge).
#' @param nfolds Cross-validation folds.
#' @param k_grid Fractional-SE tolerances to sweep in `[0, 1]`.
#' @param k Optional single tolerance; overrides the plateau search.
#' @param alpha Significance level: CI is `100(1-alpha)%`, two-sided p vs `alpha`.
#' @param intercept Fit an intercept (must be `TRUE` for binomial).
#' @param standardize glmnet column standardization. Keep `FALSE`.
#' @param seed Seed for folds; pass `NULL` inside resampling loops so the outer
#'   RNG stream is not reset.
#' @param min_vars Minimum model size (never a zero-variable model).
#' @param tailarea_tol Tolerance for the `tailarea_ok` validity flag.
#'
#' @return A list with `sweep`, `k_plateau`, `lambda`, `fit`, `table`, `n_vars`,
#'   `variables`, `tailarea`, `tailarea_ok`, `alpha`, `ci_level`.
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom selectiveInference fixedLassoInf
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom dplyr filter arrange desc
#' @export
fit_selective_lasso <- function(
    x,
    y,
    family = "binomial",
    elastic_net_alpha = 1,
    nfolds = 10,
    k_grid = seq(0, 1, 0.05),
    k = NULL,
    alpha = 0.05,
    intercept = TRUE,
    standardize = FALSE,
    seed = 141845,
    min_vars = 1,
    tailarea_tol = 0.01
) {
  
  # ---- 0. Setup + reproducible folds (NULL seed = don't reset the stream) ---
  # Seed only when asked (bootstrap passes NULL to keep the outer RNG intact)
  if (!is.null(seed)) set.seed(seed)
  # Sample size and fixed fold assignment for reproducible cross-validation
  n <- length(y)
  foldid <- sample(rep(seq_len(nfolds), length.out = n))
  
  # ---- 1. Cross-validation: lambda path + CV-error curve --------------------
  # Cross-validate the penalized logistic to get the lambda path and CV error
  cv_fit <-
    glmnet::cv.glmnet(
      x = x,
      y = y,
      family = family,
      alpha = elastic_net_alpha,
      foldid = foldid,
      standardize = standardize
    )
  
  # ---- 2. Full-path fit: single, clean source for coefficient reads ---------
  # Refit the full path standalone so exact = TRUE has a well-formed call
  gfit <-
    glmnet::glmnet(
      x = x,
      y = y,
      family = family,
      alpha = elastic_net_alpha,
      intercept = intercept,
      standardize = standardize
    )
  
  # ---- 3. Fractional-SE sweep: record the selected set at each tolerance -----
  # Index of the minimum CV error, the anchor for the SE tolerance
  imin <- which.min(cv_fit$cvm)
  # For each tolerance k, record lambda and the identity of the selected set
  sweep <-
    purrr::map_dfr(k_grid, function(k) {
      # Threshold on the CV-error scale and the most-regularized lambda under it
      thr <- cv_fit$cvm[imin] + k * cv_fit$cvsd[imin]
      lam <- max(cv_fit$lambda[cv_fit$cvm <= thr])
      # Selected (non-zero, non-intercept) variables at that lambda
      b <- as.matrix(coef(gfit, s = lam))
      vars <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
      # One row per tolerance: lambda, model size, and a set fingerprint
      tibble::tibble(
        k = k,
        lambda = lam,
        n_vars = length(vars),
        var_key = paste(sort(vars), collapse = "|")
      )
    })
  
  # ---- 4. Most stable non-empty plateau (used only when k is NULL) ----------
  # Run-length encode the fingerprints into contiguous plateaus
  runs <- rle(sweep$var_key)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  # Rank plateaus: longest first, tie-break toward more regularization
  plateaus <-
    tibble::tibble(
      var_key = runs$values,
      len = runs$lengths,
      start = starts,
      end = ends,
      n_vars = sweep$n_vars[starts]
    ) %>%
    dplyr::filter(n_vars >= min_vars) %>%
    dplyr::arrange(dplyr::desc(len), dplyr::desc(start))
  
  # ---- 5. Choose lambda: user k, else plateau center, else lambda.min -------
  # Pick lambda per priority order and record the tolerance range it spans
  if (!is.null(k)) {
    # User-supplied tolerance
    thr <- cv_fit$cvm[imin] + k * cv_fit$cvsd[imin]
    lambda <- max(cv_fit$lambda[cv_fit$cvm <= thr])
    k_plateau <- c(k, k)
  } else if (nrow(plateaus) > 0) {
    # Center of the most stable plateau
    best <- plateaus[1, ]
    mid_idx <- floor((best$start + best$end) / 2)
    lambda <- sweep$lambda[mid_idx]
    k_plateau <- c(sweep$k[best$start], sweep$k[best$end])
  } else {
    # Fallback to lambda.min when no non-empty plateau exists
    lambda <- cv_fit$lambda.min
    k_plateau <- c(NA_real_, NA_real_)
  }
  # Protection: relax lambda until at least min_vars variables are selected
  idx <- which.min(abs(cv_fit$lambda - lambda))
  while (cv_fit$nzero[idx] < min_vars && idx < length(cv_fit$lambda)) {
    idx <- idx + 1L
  }
  lambda <- cv_fit$lambda[idx]
  
  # ---- 6. Exact coefficients at the chosen lambda ---------------------------
  # Extract the exact coefficient vector at the chosen lambda (keep intercept)
  beta <- coef(object = gfit, x = x, y = y, s = lambda, exact = TRUE)
  
  # ---- 7. Selective inference (lambda * n rescales glmnet's 1/n penalty) -----
  # Selective p-values and truncated intervals at the matching standard-scale penalty
  out <-
    selectiveInference::fixedLassoInf(
      x = x,
      y = y,
      beta = beta,
      lambda = lambda * n,
      alpha = alpha,
      family = family,
      intercept = intercept
    )
  
  # ---- 8. Results table on the OR scale -------------------------------------
  # Build the OR table (two-sided + one-sided p) when anything was selected
  if (length(out$vars) > 0) {
    table <-
      tibble::tibble(
        term = names(out$vars),
        odds_ratio = exp(out$coef0),
        # out$pv is the ONE-SIDED pivot (uniform under H0). The two-sided p is
        # 2 * min(p, 1 - p), NOT 2 * p: the latter wrongly caps at 1 whenever the
        # selective estimate lands on the far side of the null (p > 0.5).
        pval = 2 * pmin(out$pv, 1 - out$pv),
        pval_1sided = out$pv,
        ci_lo = exp(out$ci[, 1]),
        ci_hi = exp(out$ci[, 2])
      ) %>%
      dplyr::arrange(pval)
    # Record the CI level so the plot can label itself
    attr(table, "ci_level") <- 1 - alpha
  } else {
    # Empty table when the model selected nothing
    table <- tibble::tibble()
  }
  # Number of selected variables
  n_vars <- length(out$vars)
  
  # ---- 9. Validity diagnostic (report only; never stop) ---------------------
  # Flag whether every achieved tail area is within tolerance of alpha/2
  if (n_vars > 0) {
    tailarea_ok <- all(abs(out$tailarea - alpha / 2) <= tailarea_tol)
  } else {
    tailarea_ok <- NA
  }
  # Print pass/fail without ever halting the run
  message(
    "Selective-inference validity check: ",
    if (isTRUE(tailarea_ok)) {
      "PASSED"
    } else if (isFALSE(tailarea_ok)) {
      "FAILED (tailarea off target - treat CIs with caution)"
    } else {
      "N/A (no variables selected)"
    }
  )
  
  # ---- 10. Return (elements in creation order) ------------------------------
  # Return the sweep, chosen lambda, raw fit, table and diagnostics
  list(
    sweep = sweep,
    k_plateau = k_plateau,
    lambda = lambda,
    fit = out,
    table = table,
    n_vars = n_vars,
    variables = names(out$vars),
    tailarea = out$tailarea,
    tailarea_ok = tailarea_ok,
    alpha = alpha,
    ci_level = 1 - alpha
  )
}


# ---- Method 2: bootstrap of the selective-inference LASSO -------------------

#' Bootstrap the selective-inference LASSO and store the raw draws
#'
#' @description
#' Resample rows with replacement `B` times; each time re-run the whole
#' select-and-infer pipeline and store the selected variables, their log-ORs,
#' and the native one-sided selective p. Summarizing is separated
#' ([summarize_bootstrap()]) so the expensive fits run once.
#'
#' @param x,y Centered design matrix and 0/1 outcome.
#' @param B Number of bootstrap replicates.
#' @param alpha Significance level passed to [fit_selective_lasso()].
#' @param k Tolerance for lambda (`0` = lambda.min, `NULL` = plateau).
#' @param elastic_net_alpha,nfolds,standardize Passed to [fit_selective_lasso()].
#' @param seed Seed for the resampling. Reproducible and invariant to the number
#'   of parallel workers, but the values differ from a plain sequential loop
#'   because parallel-safe RNG uses an independent stream per replicate.
#' @return List with `draws` (long: rep, term, log_or, si_p), `B`, `B_valid`,
#'   `n`, `alpha`, `terms`.
#' @details
#' Replicates run through [furrr::future_map()], so they parallelize across
#' whatever `future::plan()` is active (e.g. `future::plan("multisession")`);
#' the default sequential plan runs them serially.
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom furrr future_map furrr_options
#' @export
bootstrap_lasso_raw <- function(
    x,
    y,
    B = 1000,
    alpha = 0.05,
    k = 0,
    elastic_net_alpha = 1,
    nfolds = 10,
    standardize = FALSE,
    seed = 141845
) {
  # Seed the per-replicate RNG streams (reproducible, worker-count invariant)
  set.seed(seed)
  n <- length(y)
  
  # One resample-and-infer per replicate, in parallel under the active plan
  res_list <-
    furrr::future_map(
      seq_len(B),
      function(b) {
        # Draw a resample with replacement
        idx <- sample.int(n, n, replace = TRUE)
        xb <- x[idx, , drop = FALSE]
        yb <- y[idx]
        # Skip degenerate resamples with a single outcome class
        if (length(unique(yb)) < 2) {
          return(list(valid = FALSE, draw = NULL))
        }
        # Re-run the full select-and-infer pipeline (seed = NULL uses the
        # replicate's own furrr RNG stream)
        res <- tryCatch(
          suppressMessages(suppressWarnings(
            fit_selective_lasso(
              x = xb,
              y = yb,
              alpha = alpha,
              k = k,
              elastic_net_alpha = elastic_net_alpha,
              nfolds = nfolds,
              standardize = standardize,
              seed = NULL
            )
          )),
          error = function(e) NULL
        )
        # Errored fit: not a usable replicate
        if (is.null(res)) {
          return(list(valid = FALSE, draw = NULL))
        }
        # Usable but nothing selected: valid, no draw
        if (length(res$fit$vars) == 0) {
          return(list(valid = TRUE, draw = NULL))
        }
        # Store this replicate's selected variables, log-ORs and one-sided p
        draw <-
          tibble::tibble(
            rep = b,
            term = names(res$fit$vars),
            log_or = as.numeric(res$fit$coef0),
            si_p = as.numeric(res$fit$pv)
          )
        list(valid = TRUE, draw = draw)
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  
  # Collect validity flags and stack the draws into one long tibble
  valid <- vapply(res_list, function(r) r$valid, logical(1))
  draws_all <- dplyr::bind_rows(lapply(res_list, function(r) r$draw))
  # Report how many resamples were usable (stderr, so it streams under knitr)
  cat(
    "Bootstrap complete: ", sum(valid), " / ", B, " valid resamples.\n",
    sep = "", file = stderr()
  )
  # Return the draws plus book-keeping needed for summarizing
  list(
    draws = draws_all,
    B = B,
    B_valid = sum(valid),
    n = n,
    alpha = alpha,
    terms = sort(unique(draws_all$term))
  )
}


# ---- Method 3: multi-sample splitting with cross-fitting --------------------

#' Multi-sample splitting with cross-fitting (selection-corrected, on the OR scale)
#'
#' @description
#' Repeatedly split the data in half, **select** variables by LASSO on one half,
#' and **estimate + test** them by Firth-penalized logistic regression on the
#' other half. Because selection and inference never share rows, the per-split
#' p-values are valid; aggregating them over many splits gives selection-aware
#' p-values per variable on the familiar OR scale (Meinshausen, Meier & Buhlmann
#' 2009, JASA). Cross-fitting uses each half for both roles, so no data is
#' wasted. Estimation is Firth-penalized because the half-samples are small and
#' separation-prone, where ordinary MLE would blow up.
#'
#' @details
#' Splits are **stratified by the outcome**: half of each class goes to each
#' side, so both halves keep the full-sample event rate. This maximizes and
#' stabilizes the minority-class count per half — important for rare outcomes,
#' where a purely random split routinely leaves a selection half with too few
#' events for a stable LASSO fit (glmnet's "fewer than 8 observations" warning).
#'
#' Everything is built on the **cross-fit combined estimate**: within a split the
#' two directions estimate on disjoint halves, so they are independent and are
#' inverse-variance combined into one full-sample-precision estimate `beta_s` and
#' SE `se_s` for that split (two half-sample SEs combine to about the full-sample
#' SE — this is what recovers power).
#'
#' Six p-values are reported, split across two PARADIGMS. Names encode the
#' paradigm and the aggregator, so a column never implies a relationship it does
#' not have. In particular the Meinshausen family is NOT the estimate family
#' "with a penalty applied" — the two combine different things, and either may be
#' the smaller.
#'
#' **Estimate paradigm** — pool `beta_s`/`se_s` across splits, then run ONE Wald
#' test. Each block is self-contained: its OR, SE, CI and p come from the same
#' estimate, so they cohere (the CI excludes 1 exactly when p < 0.05). The two
#' blocks differ ONLY in the aggregator.
#'
#' * `OR_mean` / `SE_mean` / `ci_*_mean` / `pval_mean` — the mean split estimate;
#'   total variance = mean within-split variance + `(1 + 1/S)` times the
#'   between-split variance (a Rubin-style decomposition). Both terms are means,
#'   so it is sensitive to the splits in which only one half selected the
#'   variable: those are half-sample fits, and they drag a mean.
#' * `pval_mean_imputed` — `pval_mean`, but absent splits contribute a null
#'   effect (`beta_s = 0`), shrinking the variable toward OR = 1 in proportion to
#'   how often the LASSO dropped it.
#' * `OR_median` / `SE_median` / `ci_*_median` / `pval_median` — the median rule
#'   recommended by Chernozhukov et al. (2018) for a repeated cross-fit:
#'   `theta = median(beta_s)`, `se^2 = median(se_s^2 + (beta_s - theta)^2)`.
#'   Robust to those half-sample splits. **This is what the forest plot shows.**
#'   NOTE we borrow the aggregator, not the setting: in their work every split
#'   estimates the SAME parameter, whereas here the selected adjustment set
#'   changes from split to split, so their guarantees do not transfer. Do not
#'   describe this as DML in a manuscript.
#'   Expect `SE_median < SE_mean` as a rule, not as a finding: the per-split
#'   quantity `se_s^2 + (beta_s - theta)^2` is right-skewed, so its median sits
#'   below its mean by construction. The reassurance that the median is the right
#'   summary comes from the ORs barely moving, not from the p-values falling.
#'
#' **P-value paradigm** — compute a p in every split, then combine the p-values
#' by Meinshausen's quantile rule `min(1, quantile(p_s / gamma, gamma))`. No
#' estimate is pooled, so this family has NO OR and NO CI; there is nothing to
#' exponentiate, and nothing for an interval to agree with.
#'
#' * `pval_meinshausen` — the quantile rule over the PRESENT splits only.
#' * `pval_meinshausen_imputed` — the same over ALL splits, with `p_s = 1` for
#'   absent splits (Meinshausen's own convention).
#' * `pval_meinshausen_fwer` — faithful Meinshausen, Meier & Buhlmann (2009):
#'   additionally Bonferroni each split's p by that split's selected-model size
#'   `|S^(b)|` before aggregating. This alone is **FWER-adjusted** across all
#'   candidate predictors, so it is NOT comparable with the unadjusted p-values
#'   of every other method — compare it to alpha directly. (We use a fixed
#'   `gamma`; MMB also offer an adaptive search over `gamma` with a further
#'   `(1 - log gamma_min)` factor.)
#'
#' The Meinshausen p-values are typically SMALLER than the estimate-paradigm
#' ones, because a per-split p carries only that split's within-variance — the
#' between-split spread never enters. That is not a penalty going missing; it is
#' the two paradigms answering different questions.
#'
#' Neither pooled OR/CI has a formal coverage proof, because different splits may
#' select different adjustment sets; read them as the typical effect across the
#' selected models.
#'
#' @param x,y Centered design matrix and 0/1 outcome.
#' @param n_splits Number of stratified half-splits (each contributes two
#'   cross-fitted directions).
#' @param gamma Quantile for p-value aggregation (`0.5` = median rule).
#' @param nfolds CV folds for the LASSO selection step.
#' @param elastic_net_alpha glmnet mixing (1 = lasso).
#' @param standardize glmnet column standardization. Keep `FALSE`.
#' @param seed Seed for the splits. Reproducible and invariant to the number of
#'   parallel workers, but the values differ from a plain sequential loop
#'   because parallel-safe RNG uses an independent stream per split.
#'
#' @return Tibble of `term`, `n_chosen`, `pct_chosen` (selection frequency); the
#'   mean-aggregated estimate (`OR_mean`, `SE_mean`, `ci_lo_mean`, `ci_hi_mean`,
#'   `pval_mean`, `pval_mean_imputed`); the median-aggregated estimate
#'   (`OR_median`, `SE_median`, `ci_lo_median`, `ci_hi_median`, `pval_median`);
#'   and the three Meinshausen p-values (`pval_meinshausen`,
#'   `pval_meinshausen_imputed`, `pval_meinshausen_fwer`), which have no OR or CI
#'   by construction. All are `NA` for variables present in fewer than 2 splits.
#' @details
#' Splits run through [furrr::future_map()], so they parallelize across whatever
#' `future::plan()` is active; the default sequential plan runs them serially.
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
#' @importFrom dplyr arrange desc
#' @importFrom furrr future_map furrr_options
#' @export
fit_split_lasso <- function(
    x,
    y,
    n_splits = 1000,
    gamma = 0.5,
    nfolds = 10,
    elastic_net_alpha = 1,
    standardize = FALSE,
    seed = 141845
) {
  # Seed the per-split RNG streams (reproducible, worker-count invariant)
  set.seed(seed)
  n <- length(y)
  terms <- colnames(x)
  
  # Storage: two cross-fitted directions per split. We keep the held-out log-OR
  # and its SE per direction; every p-value is derived later from the cross-fit
  # COMBINED estimate, so no per-direction p-values are stored.
  n_dir <- 2L * n_splits
  bmat <- matrix(
    NA_real_, nrow = n_dir, ncol = length(terms), dimnames = list(NULL, terms)
  )
  semat <- matrix(
    NA_real_, nrow = n_dir, ncol = length(terms), dimnames = list(NULL, terms)
  )
  sel_mat <- matrix(
    FALSE, nrow = n_dir, ncol = length(terms), dimnames = list(NULL, terms)
  )
  
  # Row indices grouped by outcome class, reused for every stratified split
  idx_by_class <- split(seq_len(n), y)
  
  # One stratified split per iteration, in parallel under the active plan. Each
  # returns its two cross-fitted directions (NULL where a direction was
  # single-class, selected nothing, or the held-out fit failed).
  split_results <-
    furrr::future_map(
      seq_len(n_splits),
      function(s) {
        # Stratified half/half partition: draw half of EACH class into part_a so
        # both halves keep the full-sample event rate (maximizes and stabilizes
        # the minority-class count per half).
        part_a <- unlist(lapply(idx_by_class, function(ix) {
          sample(ix, floor(length(ix) / 2))
        }))
        part_b <- setdiff(seq_len(n), part_a)
        # Cross-fit: (select A, test B) then (select B, test A)
        dirs <- list(
          list(sel = part_a, est = part_b),
          list(sel = part_b, est = part_a)
        )
        lapply(dirs, function(dir) {
          y_sel <- y[dir$sel]
          y_est <- y[dir$est]
          # Skip if either half is single-class (can't select or test)
          if (length(unique(y_sel)) < 2 || length(unique(y_est)) < 2) {
            return(NULL)
          }
          # Select by LASSO on the selection half (glmnet's "fewer than 8
          # observations" warning is expected here, so suppress it)
          vars <- tryCatch(
            suppressWarnings(
              select_lasso_vars(
                x = x[dir$sel, , drop = FALSE],
                y = y_sel,
                nfolds = nfolds,
                elastic_net_alpha = elastic_net_alpha,
                standardize = standardize
              )
            ),
            error = function(e) character()
          )
          # Nothing selected -> contributes only the default p = 1
          if (length(vars) == 0) return(NULL)
          # Estimate + test on the held-out half, Firth-penalized so the small,
          # separation-prone half stays finite
          est <- firth_estimate(x[dir$est, , drop = FALSE], y_est, vars)
          if (is.null(est)) return(NULL)
          list(terms = est$terms, coef = est$coef, se = est$se)
        })
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  
  # Flatten to per-direction results (2 * n_splits, NULLs kept) and fill the
  # matrices; row order is irrelevant to the order-free aggregation below
  dir_results <- do.call(c, split_results)
  for (d in seq_along(dir_results)) {
    r <- dir_results[[d]]
    if (is.null(r)) next
    bmat[d, r$terms] <- r$coef
    semat[d, r$terms] <- r$se
    sel_mat[d, r$terms] <- TRUE
  }
  
  # Per-split selected-model size: the DISTINCT variables tested in that split
  # (union of its two cross-fit directions). This is Meinshausen's |S^(b)|, the
  # Bonferroni factor used by `pval_meinshausen_fwer` below.
  msize <- vapply(split_results, function(dd) {
    tt <- unlist(lapply(dd, function(r) if (is.null(r)) character() else r$terms))
    length(unique(tt))
  }, integer(1))
  
  # Meinshausen quantile rule at the chosen gamma (gamma = 0.5 -> 2 * median)
  meinshausen <- function(p) {
    min(1, quantile(p / gamma, probs = gamma, names = FALSE))
  }
  
  # Aggregate across splits into one row per variable
  purrr::map_dfr(terms, function(v) {
    n_chosen <- sum(sel_mat[, v])
    pct_chosen <- round(100 * n_chosen / n_dir, 1)
    # Cross-fit combine. Within a split the two directions estimate on DISJOINT
    # halves, so they are independent and are inverse-variance combined into one
    # full-sample-precision estimate. Reshape to 2 x n_splits (column = split).
    bm <- matrix(bmat[, v], nrow = 2)
    sm <- matrix(semat[, v], nrow = 2)
    km <- matrix(sel_mat[, v], nrow = 2)
    w <- ifelse(km, 1 / sm^2, 0)
    w[!is.finite(w)] <- 0
    b0 <- ifelse(km, bm, 0)
    b0[!is.finite(b0)] <- 0
    wsum <- colSums(w)
    present <- wsum > 0       # splits where v entered at least one direction
    S <- sum(present)
    # Too few splits to pool anything meaningful
    if (S < 2) {
      return(tibble::tibble(
        term = v, n_chosen = n_chosen, pct_chosen = pct_chosen,
        OR_mean = NA_real_, SE_mean = NA_real_,
        ci_lo_mean = NA_real_, ci_hi_mean = NA_real_,
        pval_mean = NA_real_, pval_mean_imputed = NA_real_,
        OR_median = NA_real_, SE_median = NA_real_,
        ci_lo_median = NA_real_, ci_hi_median = NA_real_,
        pval_median = NA_real_,
        pval_meinshausen = NA_real_, pval_meinshausen_imputed = NA_real_,
        pval_meinshausen_fwer = NA_real_
      ))
    }
    beta_s <- colSums(w * b0)[present] / wsum[present]
    se_s <- 1 / sqrt(wsum[present])
    
    # ---- Estimate paradigm: pool the split estimates, then one Wald test -----
    # Both blocks below pool `beta_s`/`se_s` across splits and run ONE Wald test.
    # Each is self-contained: its OR, SE, CI and p come from the same estimate,
    # so they cohere. The two differ only in the aggregator, mean vs median.
    #
    # mean: the REPORTED effect and the one in the forest plot. Total variance =
    # mean within-split variance + (1 + 1/S) times the between-split variance
    # (a Rubin-style decomposition; both terms are means, hence outlier-prone).
    est_mean <- mean(beta_s)
    se_mean <- sqrt(mean(se_s^2) + (1 + 1 / S) * stats::var(beta_s))
    pval_mean <- 2 * stats::pnorm(-abs(est_mean / se_mean))
    or_mean <- exp(est_mean)
    ci_lo_mean <- exp(est_mean - stats::qnorm(0.975) * se_mean)
    ci_hi_mean <- exp(est_mean + stats::qnorm(0.975) * se_mean)
    # mean_imputed: absent splits contribute a null effect (log-OR = 0) with the
    # precision a present split typically had, shrinking v toward OR = 1 in
    # proportion to how often the LASSO dropped it
    b_all <- rep(0, n_splits)
    b_all[present] <- beta_s
    se_all <- rep(mean(se_s), n_splits)
    se_all[present] <- se_s
    est_i <- mean(b_all)
    se_i <- sqrt(mean(se_all^2) + (1 + 1 / n_splits) * stats::var(b_all))
    pval_mean_imputed <- 2 * stats::pnorm(-abs(est_i / se_i))
    # median: the aggregator recommended by Chernozhukov et al. (2018) for a
    # repeated cross-fit. NOTE we borrow the aggregator, not the setting: in
    # their work every split estimates the SAME parameter, whereas here the
    # selected adjustment set changes from split to split, so their guarantees
    # do not transfer. Two differences from `mean`: the median split estimate,
    # which the half-sample splits (only one half selected v) cannot drag; and
    # the between-split spread entering INSIDE the median, per split, rather
    # than as a separate Rubin term.
    est_median <- stats::median(beta_s)
    se_median <- sqrt(stats::median(se_s^2 + (beta_s - est_median)^2))
    pval_median <- 2 * stats::pnorm(-abs(est_median / se_median))
    or_median <- exp(est_median)
    ci_lo_median <- exp(est_median - stats::qnorm(0.975) * se_median)
    ci_hi_median <- exp(est_median + stats::qnorm(0.975) * se_median)
    
    # ---- p-value paradigm: combine the per-split p-values (Meinshausen) ------
    # No estimate is pooled here, so this family has no OR and no CI - there is
    # nothing to exponentiate. That is why it cannot be lined up against a CI.
    # Each split's own Wald p, against that split's own SE
    p_s <- 2 * stats::pnorm(-abs(beta_s / se_s))
    # meinshausen: quantile rule over the present splits only
    pval_meinshausen <- meinshausen(p_s)
    # meinshausen_imputed: absent splits contribute p = 1 (Meinshausen's rule)
    p_all <- rep(1, n_splits)
    p_all[present] <- p_s
    pval_meinshausen_imputed <- meinshausen(p_all)
    # meinshausen_fwer: faithful Meinshausen, Meier & Buhlmann - additionally
    # Bonferroni each split's p by that split's selected-model size, so the
    # result is FWER-adjusted across variables. Compare it to alpha directly;
    # it is NOT comparable to the unadjusted p-values of the other methods.
    p_fwer <- rep(1, n_splits)
    p_fwer[present] <- pmin(1, p_s * msize[present])
    pval_meinshausen_fwer <- meinshausen(p_fwer)
    
    tibble::tibble(
      term = v,
      n_chosen = n_chosen,
      pct_chosen = pct_chosen,
      OR_mean = or_mean,
      SE_mean = se_mean,
      ci_lo_mean = ci_lo_mean,
      ci_hi_mean = ci_hi_mean,
      pval_mean = pval_mean,
      pval_mean_imputed = pval_mean_imputed,
      OR_median = or_median,
      SE_median = se_median,
      ci_lo_median = ci_lo_median,
      ci_hi_median = ci_hi_median,
      pval_median = pval_median,
      pval_meinshausen = pval_meinshausen,
      pval_meinshausen_imputed = pval_meinshausen_imputed,
      pval_meinshausen_fwer = pval_meinshausen_fwer
    )
  }) %>%
    dplyr::arrange(dplyr::desc(pct_chosen))
}


# ---- Method 4: Firth penalized logistic -------------------------------------

#' Firth penalized logistic on a fixed predictor set (separation-robust)
#'
#' @description
#' Fit `logistf::logistf()` on a given set of variables. Firth's penalty keeps
#' estimates finite under separation and improves small-sample behavior. NOTE:
#' this corrects the *estimation*, not the *selection* — if `vars` came from a
#' data-driven search, the p-values remain post-selection (naive).
#'
#' @param x,y Design matrix and 0/1 outcome.
#' @param vars Variable names to fit.
#' @return Tibble of `term`, `OR`, `SE`, `ci`, `pval`; empty if `vars` is empty.
#' @importFrom tibble tibble
#' @export
fit_firth <- function(x, y, vars) {
  # Nothing to fit if no variables were passed
  if (length(vars) == 0) return(tibble::tibble(term = character()))
  # Build the model frame from the selected columns and the outcome
  df <- data.frame(.y = y, x[, vars, drop = FALSE], check.names = FALSE)
  # Fit the Firth-penalized logistic (returns NULL on failure)
  fit <- tryCatch(logistf::logistf(.y ~ ., data = df), error = function(e) NULL)
  # Fall back to a bare term list if the fit failed
  if (is.null(fit)) return(tibble::tibble(term = vars))
  # Indices of the non-intercept coefficients
  idx <- which(names(coef(fit)) != "(Intercept)")
  # Return ORs, SEs, profile-penalized CIs and penalized LR p-values
  tibble::tibble(
    term = names(coef(fit))[idx],
    OR = exp(coef(fit)[idx]),
    SE = sqrt(diag(vcov(fit)))[idx],
    ci = ci_str(exp(fit$ci.lower[idx]), exp(fit$ci.upper[idx])),
    pval = fit$prob[idx]
  )
}


# ---- Method 5: forward stepwise (native fs + logistic refit) ----------------

#' Forward stepwise selective inference, plus a logistic refit for ORs
#'
#' @description
#' Run forward stepwise with AIC stopping and selective inference on the LINEAR
#' scale via [selectiveInference::fsInf()] (the package supports Gaussian only),
#' then refit an ordinary logistic model on the selected variables to obtain
#' interpretable odds ratios. Two views are returned: the native (valid, linear)
#' inference and the logistic ORs (post-selection, naive).
#'
#' @param x,y Centered design matrix and 0/1 outcome (treated as continuous by fs).
#' @param maxsteps Maximum forward-selection steps.
#' @param seed Seed for the fs fit.
#' @return List with `fsinf` (raw), `native` (linear coef/CI/selective p) and
#'   `logistic` (OR/SE/CI/p) tibbles.
#' @importFrom selectiveInference fs fsInf estimateSigma
#' @importFrom tibble tibble
#' @export
fit_stepwise <- function(x, y, maxsteps = 2000, seed = 141845) {
  # Seed the forward-selection path
  set.seed(seed)
  # Run forward stepwise (linear/Gaussian) without internal normalization
  fsfit <-
    selectiveInference::fs(
      x = x,
      y = y,
      maxsteps = maxsteps,
      intercept = TRUE,
      normalize = FALSE
    )
  # Estimate residual sigma for the selective inference
  sigmahat <-
    selectiveInference::estimateSigma(
      x = x,
      y = y,
      intercept = TRUE,
      standardize = FALSE
    )$sigmahat
  # Selective inference at the AIC stopping point
  fsinf <-
    selectiveInference::fsInf(fsfit, type = "aic", sigma = sigmahat)
  
  # Map the selected variable indices to names
  vars <- colnames(x)[fsinf$vars]
  
  # Native (linear-scale) view: coefficient, selective CI and selective p.
  # fsInf's pv is the ONE-SIDED pivot (uniform under H0), so the two-sided p is
  # 2 * min(p, 1 - p) — same convention as fit_selective_lasso().
  native <-
    tibble::tibble(
      term = vars,
      coef = as.numeric(fsinf$vmat %*% fsinf$y),
      ci = ci_str(fsinf$ci[, 1], fsinf$ci[, 2]),
      pval = 2 * pmin(fsinf$pv, 1 - fsinf$pv),
      pval_1sided = fsinf$pv
    )
  
  # Logistic refit on the selected variables for interpretable ORs
  est <- glm_estimate(x, y, vars)
  # Wald CIs for the same logistic refit (NULL if the estimate failed)
  ci <- if (!is.null(est)) stats::confint.default(glm(
    .y ~ ., data = data.frame(.y = y, x[, vars, drop = FALSE], check.names = FALSE),
    family = binomial()
  )) else NULL
  # Build the logistic OR view (or a bare term list on failure)
  logistic <-
    if (is.null(est)) {
      tibble::tibble(term = vars)
    } else {
      tibble::tibble(
        term = est$terms,
        OR = exp(est$coef),
        SE = est$se,
        ci = ci_str(exp(ci[est$terms, 1]), exp(ci[est$terms, 2])),
        pval = est$pval
      )
    }
  
  # Return the raw object and both views
  list(fsinf = fsinf, native = native, logistic = logistic)
}


# ---- Method 6: univariable logistic -----------------------------------------

#' Univariable logistic regression across many predictors
#'
#' @description
#' Fit one logistic model per predictor (optionally adjusted for a fixed set),
#' with robust (sandwich) CIs, and stack the results.
#'
#' @param df Data frame.
#' @param outcome Outcome column name.
#' @param cols Predictor column names to test one at a time.
#' @param adjust Character vector of always-included adjustors (or `NULL`).
#' @param ci_high_cap Drop rows whose upper CI exceeds this (unstable fits).
#' @return Tibble in the `Predictors / Odds Ratios (OR) / SE / ...` schema.
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter arrange
#' @export
fit_univariable <- function(
    df,
    outcome,
    cols,
    adjust = NULL,
    ci_high_cap = 50
) {
  # Adjustor prefix for the formula right-hand side (empty when unadjusted)
  rhs_base <- if (is.null(adjust)) "" else paste(paste(adjust, collapse = " + "), "+")
  # Fit one robust logistic model per predictor and stack the tidy rows
  out <-
    purrr::map_dfr(cols, function(v) {
      # Skip a predictor that is itself an adjustor
      if (v %in% adjust) return(tibble::tibble())
      # Build and fit the (optionally adjusted) univariable model
      model <- as.formula(paste(outcome, "~", rhs_base, v))
      fit <- suppressWarnings(glm(model, data = df, family = binomial()))
      # Robust (sandwich) covariance and exponentiated, stylised coefficients
      robust <- sandwich::vcovHC(fit, type = "HC0")
      lmtest::coeftest(fit, vcov. = robust) %>%
        broom::tidy(conf.int = TRUE) %>%
        dplyr::mutate(dplyr::across(c(estimate, conf.low, conf.high), exp)) %>%
        dplyr::transmute(
          Predictors = term,
          `Odds Ratios (OR)` = round_or(estimate),
          SE = round(std.error, 2),
          `Z-scores` = statistic,
          `P-values` = round(p.value, 3),
          `CI (low)` = round_or(conf.low),
          `CI (high)` = round_or(conf.high)
        )
    })
  # Drop the intercept, the adjustor rows and unstable (huge-CI) rows
  out %>%
    dplyr::filter(Predictors != "(Intercept)") %>%
    dplyr::filter(
      is.null(adjust) | !grepl(paste(adjust, collapse = "|"), Predictors)
    ) %>%
    dplyr::filter(`CI (high)` < ci_high_cap) %>%
    dplyr::arrange(`P-values`)
}


# ---- Method 7: p-value-screened multivariable logistic ----------------------

#' Multivariable logistic from p-value-screened predictors, with robust CIs
#'
#' @description
#' Screen predictors from a univariable table at `pval_threshold`, drop
#' `unwanted`, fit one logistic model, and return sandwich-robust ORs. NOTE: the
#' p-values are post-selection (screened on the same data), so treat them as
#' exploratory.
#'
#' @param df Data frame.
#' @param outcome Outcome column name.
#' @param univ_tables List of univariable tibbles to screen from.
#' @param unwanted Character vector of predictors to exclude.
#' @param pval_threshold Screening threshold.
#' @param drop_na_cols Columns that must be complete before fitting.
#' @return List with `fit` (glm) and `table` (OR schema tibble).
#' @importFrom dplyr bind_rows filter pull transmute
#' @importFrom stringr str_replace
#' @export
fit_multivariable <- function(
    df,
    outcome,
    univ_tables,
    unwanted = character(),
    pval_threshold = 0.2,
    drop_na_cols = NULL
) {
  # Bind the univariable screen tables and check the columns we screen on exist
  screened <- dplyr::bind_rows(univ_tables)
  needed <- setdiff(c("Predictors", "P-values"), colnames(screened))
  if (length(needed) > 0) {
    stop(
      "fit_multivariable: screen tables lack column(s) ",
      paste(needed, collapse = ", "), ". Found: ",
      if (ncol(screened) == 0) "<none>" else paste(colnames(screened), collapse = ", "),
      call. = FALSE
    )
  }
  # Screen predictors below the threshold and clean off factor-level suffixes
  predictors <-
    screened %>%
    dplyr::filter(`P-values` < pval_threshold) %>%
    dplyr::pull(`Predictors`) %>%
    stringr::str_replace("^([a-z0-9_]*)([A-Z].*)$", "\\1") %>%
    unique()
  # Remove unwanted predictors
  predictors <- predictors[!predictors %in% unwanted]
  
  # Restrict to complete cases on the required columns when requested
  if (!is.null(drop_na_cols)) df <- tidyr::drop_na(df, dplyr::all_of(drop_na_cols))
  # Build and fit the single multivariable model
  model <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  fit <- glm(model, data = df, family = binomial())
  
  # Robust (sandwich) covariance for the reported CIs
  robust <- sandwich::vcovHC(fit, type = "HC0")
  # Style the exponentiated coefficients into the shared OR schema
  table <-
    lmtest::coeftest(fit, vcov. = robust) %>%
    broom::tidy(conf.int = TRUE) %>%
    dplyr::mutate(dplyr::across(c(estimate, conf.low, conf.high), exp)) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::transmute(
      Predictors = term,
      `Odds Ratios (OR)` = round_or(estimate),
      SE = round(std.error, 2),
      `Z-scores` = statistic,
      `P-values` = round(p.value, 3),
      `CI (low)` = round_or(conf.low),
      `CI (high)` = round_or(conf.high)
    )
  # Return the fit, the stylised table and the chosen predictors
  list(fit = fit, table = table, predictors = predictors)
}


# ---- Bootstrap summary + comparison blocks ----------------------------------

#' Summarize a raw bootstrap object (cheap; re-run freely)
#'
#' @param raw Output of [bootstrap_lasso_raw()].
#' @param alpha Percentile-CI level.
#' @return Tibble with conditional/unconditional OR, CIs, p-values, and the
#'   distribution of the native one-sided selective p.
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter arrange desc
#' @export
summarize_bootstrap <- function(raw, alpha = 0.05) {
  # Denominator and percentile-CI quantiles
  B <- raw$B_valid
  qlo <- alpha / 2
  qhi <- 1 - alpha / 2
  # Two-sided percentile p, consistent with the percentile CI (no tie-split)
  perc_p <- function(v) pmin(1, 2 * min(mean(v <= 0), mean(v >= 0)))
  # One summarized row per variable
  purrr::map_dfr(raw$terms, function(v) {
    # Draws for this variable; conditional (selected only) vs zero-padded
    d <- dplyr::filter(raw$draws, term == v)
    lo <- d$log_or
    nch <- length(lo)
    loall <- c(lo, rep(0, B - nch))
    # Conditional, unconditional and selective-p percentile summaries
    cc <- exp(quantile(lo, c(qlo, qhi), names = FALSE))
    uu <- exp(quantile(loall, c(qlo, qhi), names = FALSE))
    sp <- quantile(d$si_p, c(qlo, qhi), names = FALSE, na.rm = TRUE)
    # Assemble the row: selection frequency, ORs, CIs and p-values
    tibble::tibble(
      term = v,
      n_chosen = nch,
      pct_chosen = 100 * nch / B,
      OR_cond = exp(median(lo)),
      OR_uncond = exp(median(loall)),
      SE = sd(lo),
      ci_cond = ci_str(cc[1], cc[2]),
      ci_uncond = ci_str(uu[1], uu[2]),
      si_pval_median = median(d$si_p, na.rm = TRUE),
      si_pval_int = ci_str(sp[1], sp[2]),
      cond_pval = perc_p(lo),
      uncond_pval = perc_p(loall)
    )
  }) %>%
    dplyr::arrange(dplyr::desc(pct_chosen))
}


# ---- Combined-table adapters (one row per term, prefixed columns) -----------

#' Selective-inference block + plot table: primary + wider CI, one-sided p
#'
#' Returns `list(block, plot)`: `block` is the string-CI row for the combined
#' table; `plot` is the numeric primary table (with its `ci_level` attribute)
#' ready for [plot_forest()].
#' @keywords internal
si_block <- function(x, y, k, alpha = 0.05, nfolds = 10, seed = 141845) {
  # Wider companion interval at twice the tail (alpha 0.05 -> 95% and 90%)
  alpha_wide <- min(2 * alpha, 0.5)
  # Same seed on both fits so they share folds and differ only in CI level
  f_primary <- suppressMessages(
    fit_selective_lasso(x, y, alpha = alpha, k = k, nfolds = nfolds, seed = seed)
  )
  f_wide <- suppressMessages(
    fit_selective_lasso(x, y, alpha = alpha_wide, k = k, nfolds = nfolds, seed = seed)
  )
  # Level-labeled CI column names (e.g. "ci95" and "ci90")
  col_primary <- paste0("ci", round(100 * (1 - alpha)))
  col_wide <- paste0("ci", round(100 * (1 - alpha_wide)))
  # String-CI block for the combined table (one row per term, both CI levels)
  block <-
    f_primary$table %>%
    dplyr::transmute(
      term,
      OR = odds_ratio,
      SE = NA_real_,
      !!col_primary := ci_str(ci_lo, ci_hi),
      pval = pval,
      pval_1sided = pval_1sided
    ) %>%
    dplyr::left_join(
      f_wide$table %>% dplyr::transmute(term, !!col_wide := ci_str(ci_lo, ci_hi)),
      by = "term"
    )
  # Numeric primary table (keeps the ci_level attribute) for the forest plot
  list(block = block, plot = f_primary$table)
}

#' Bootstrap block for the combined table (internal)
#' @keywords internal
boot_block <- function(raw, alpha = 0.05) {
  # Denominator, quantiles and the percentile p helper
  B <- raw$B_valid
  qlo <- alpha / 2
  qhi <- 1 - alpha / 2
  perc_p <- function(v) pmin(1, 2 * min(mean(v <= 0), mean(v >= 0)))
  # One row per variable: n_chosen/pct first, si_pval before cond/uncond
  purrr::map_dfr(raw$terms, function(v) {
    # Conditional draws and the zero-padded unconditional version
    d <- dplyr::filter(raw$draws, term == v)
    lo <- d$log_or
    nch <- length(lo)
    loall <- c(lo, rep(0, B - nch))
    # Percentile CIs (conditional, unconditional, selective p)
    cc <- exp(quantile(lo, c(qlo, qhi), names = FALSE))
    uu <- exp(quantile(loall, c(qlo, qhi), names = FALSE))
    sp <- quantile(d$si_p, c(qlo, qhi), names = FALSE, na.rm = TRUE)
    # Assemble the block row
    tibble::tibble(
      term = v,
      n_chosen = nch,
      pct_chosen = 100 * nch / B,
      OR_cond = exp(median(lo)),
      OR_uncond = exp(median(loall)),
      SE = sd(lo),
      ci_cond = ci_str(cc[1], cc[2]),
      ci_uncond = ci_str(uu[1], uu[2]),
      si_pval_median = median(d$si_p, na.rm = TRUE),
      si_pval_int = ci_str(sp[1], sp[2]),
      cond_pval = perc_p(lo),
      uncond_pval = perc_p(loall)
    )
  })
}

#' Split block + plot table for the combined table (internal)
#'
#' Returns `list(block, plot)`: `block` folds the CI bounds into one string
#' column for the combined table; `plot` is the numeric table (95% percentile
#' CI) in the shared forest schema, ready for [plot_forest()].
#' @keywords internal
split_block <- function(
    x,
    y,
    n_splits = 1000,
    nfolds = 10,
    seed = 141845,
    min_pct_chosen = 5
) {
  # Run the split method once; derive both the block and the plot table from it
  raw <- fit_split_lasso(x, y, n_splits = n_splits, nfolds = nfolds, seed = seed)
  # Block for the CSV: EVERY variable, however rarely selected (no filtering)
  block <-
    raw %>%
    dplyr::transmute(
      term,
      n_chosen = n_chosen,
      pct_chosen = pct_chosen,
      # Estimate paradigm, mean aggregator. Self-contained: OR, SE, CI and p
      # all come from the one estimate. This is what the forest plot shows.
      OR_mean = OR_mean,
      SE_mean = SE_mean,
      ci_mean = ci_str(ci_lo_mean, ci_hi_mean),
      pval_mean = pval_mean,
      pval_mean_imputed = pval_mean_imputed,
      # Estimate paradigm, median aggregator. Also self-contained, and never
      # mixed with the mean columns above.
      OR_median = OR_median,
      SE_median = SE_median,
      ci_median = ci_str(ci_lo_median, ci_hi_median),
      pval_median = pval_median,
      # P-value paradigm. No estimate is pooled, so there is no OR and no CI to
      # pair these with. `_fwer` is the only FWER-adjusted p in the table.
      pval_meinshausen = pval_meinshausen,
      pval_meinshausen_imputed = pval_meinshausen_imputed,
      pval_meinshausen_fwer = pval_meinshausen_fwer
    )
  # Plot table: drop variables selected in fewer than min_pct_chosen % of
  # directions - their intervals are noise and would swamp the forest.
  # The plot shows the MEDIAN aggregation (Chernozhukov et al. 2018): the median
  # is robust to the splits in which only one half selected the variable, which
  # are fitted on half the patients and drag a mean. OR, CI and p are all taken
  # from that one estimate, so the three cohere.
  plot <-
    raw %>%
    dplyr::filter(pct_chosen >= min_pct_chosen) %>%
    dplyr::transmute(
      term,
      odds_ratio = OR_median,
      pval = pval_median,
      ci_lo = ci_lo_median,
      ci_hi = ci_hi_median
    )
  attr(plot, "ci_level") <- 0.95
  list(block = block, plot = plot)
}


# ---- The compute orchestrator (computes only, never writes) -----------------

#' Compute every enabled method across all conditions (no file writing)
#'
#' @description
#' Runs each enabled method end-to-end (bootstrap and stepwise included) and
#' returns a named list of tidy tables: `combined` (the wide table, one row per
#' condition/term with all methods side by side), one un-prefixed table per
#' method, and `raw_all` (the bootstrap draws, kept so you can save them). It
#' writes NOTHING — pass the result to [write_all_methods()]. Everything is
#' driven by `methods`, so a single call is the only entry point.
#'
#' Univariable and multivariable are the least-interesting models, so their
#' inputs are bundled to keep this signature lean. Univariable is a standalone
#' computation done elsewhere and passed in as one table (`uni_table`), attached
#' to the `all` rows. Multivariable IS computed here, once per condition: it
#' screens the univariable predictors (`multi_params$screen_tables`) at
#' `multi_params$pval_threshold`, drops that condition's `multi_unwanted` set,
#' and fits one robust logistic model. A condition with `multi_unwanted = NULL`
#' gets no multivariable (e.g. the collinear `all` and the `_no_surgery` designs).
#'
#' @param conditions Named list; each element has `x`, `y` (the centered design
#'   every per-condition LASSO method runs across) and, when a multivariable is
#'   wanted, `multi_unwanted` (columns that condition drops from the screen;
#'   `NULL` skips the multivariable for that condition).
#' @param uni_table Precomputed univariable tibble (stylized schema), attached
#'   once to the `all` rows; `NULL` to skip.
#' @param multi_params List of shared multivariable settings: `df`, `outcome`,
#'   `screen_tables` (univariable tables to screen predictors from),
#'   `pval_threshold` (default 0.2) and `drop_na_cols`. `NULL` skips it.
#' @param label Name of the outcome being fitted, printed in the opening status
#'   banner so a multi-outcome batch is readable. Defaults to
#'   `multi_params$outcome`.
#' @param methods Named logical list toggling which methods to run; missing
#'   entries default to `TRUE`.
#' @param si_alpha Selective-inference significance level (95% CI at 0.05); the
#'   companion wider CI is reported at `2 * si_alpha` (90% at 0.05).
#' @param seed Master RNG seed. It is RE-APPLIED (`set.seed`) at the start of
#'   every stochastic method and every condition, so each method is reproducible
#'   and independent of the order methods run in (i.e. the seed IS reset per
#'   method, not set once for the whole run). The lone exception is the
#'   per-replicate fit inside the bootstrap, which runs with `seed = NULL` so the
#'   resampling loop's own stream is not reset mid-way.
#' @param nfolds,n_splits CV folds and number of random splits.
#' @param split_min_pct Minimum `pct_chosen` for a variable to appear in the
#'   split FOREST PLOTS; the split CSV always keeps every variable.
#' @param boot_B,boot_alpha,boot_k Bootstrap replicates, CI level and lambda rule.
#' @param verbose Print an opening banner naming the outcome, then a timestamped
#'   "<condition> / <method>: started/finished" line around each fit. Written to
#'   stderr, not via `message()`, so it streams to the console under
#'   `rmarkdown::render()` instead of being captured into the document. Set
#'   `FALSE` to silence.
#' @return Named list of tibbles (`combined` + one per enabled method), plus
#'   `raw_all` (bootstrap draws) and `plot_tables` (numeric forest tables per
#'   condition/method) — neither of the last two is written as a CSV.
#' @importFrom purrr imap imap_dfr reduce map
#' @importFrom dplyr bind_rows full_join mutate arrange desc
#' @export
compute_all_methods <- function(
    conditions,
    uni_table = NULL,
    multi_params = NULL,
    label = NULL,
    methods = NULL,
    si_alpha = 0.05,
    seed = 141845,
    nfolds = 10,
    n_splits = 1000,
    split_min_pct = 5,
    boot_B = 1000,
    boot_alpha = 0.05,
    boot_k = 0,
    verbose = TRUE
) {
  
  # ---- Resolve toggles (any missing entry defaults to TRUE) -----------------
  # Merge the caller's toggles over an all-TRUE default
  m <- utils::modifyList(
    list(
      selective_lasso = TRUE,
      bootstrap = TRUE,
      split = TRUE,
      firth = TRUE,
      stepwise = TRUE,
      univariable = TRUE,
      multivariable = TRUE
    ),
    if (is.null(methods)) list() else methods
  )
  
  # ---- Status helpers (timestamped, condition-aware; verbose toggles them) --
  # Print one "<condition> / <method>: <state>" line. We cat() straight to
  # stderr rather than using message(): knitr CAPTURES message() into the
  # rendered document, so during rmarkdown::render() it never reaches the
  # console (and is dropped entirely if the chunk sets message = FALSE).
  # Writing to stderr bypasses that capture and streams live in both modes.
  announce <- function(cond, method, state) {
    if (verbose) {
      cat(
        format(Sys.time(), "%H:%M:%S"), "  [", cond, "] ", method, ": ", state,
        "\n",
        sep = "", file = stderr()
      )
      flush(stderr())
    }
  }
  # Every method that fails, as "[condition] method", reported at the end
  failures <- character()
  # Wrap a step with a started / finished (+ elapsed seconds) status pair.
  # A method that ERRORS is not fatal: the failure is recorded, NULL is returned,
  # and the run carries on. One method that cannot fit one condition (forward
  # stepwise on a near-collinear design, say) must not destroy every other
  # method's results for that outcome — results that can take minutes to compute.
  # Downstream, `b$<method> <- NULL` simply adds no element, so a failed method
  # drops out of the combined table rather than corrupting it.
  timed <- function(cond, method, thunk) {
    announce(cond, method, "started")
    t0 <- Sys.time()
    out <- tryCatch(
      thunk(),
      error = function(e) {
        announce(cond, method, paste0("FAILED — ", conditionMessage(e)))
        failures <<- c(failures, paste0("[", cond, "] ", method))
        NULL
      }
    )
    secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    if (!is.null(out)) {
      announce(cond, method, paste0("finished (", secs, "s)"))
    }
    out
  }
  
  # ---- Announce the run ------------------------------------------------------
  # Name the outcome up front: in a multi-outcome batch the per-method status
  # lines below are meaningless without it. Falls back to the multivariable
  # outcome column when no explicit label is passed.
  run_label <- if (!is.null(label)) {
    label
  } else if (!is.null(multi_params$outcome)) {
    multi_params$outcome
  } else {
    "<unnamed outcome>"
  }
  if (verbose) {
    cat(
      "\n", strrep("=", 72), "\n",
      format(Sys.time(), "%H:%M:%S"), "  FITTING OUTCOME: ", run_label, "\n",
      "  conditions: ", paste(names(conditions), collapse = ", "), "\n",
      "  methods:    ", paste(names(m)[unlist(m)], collapse = ", "), "\n",
      strrep("=", 72), "\n",
      sep = "", file = stderr()
    )
    flush(stderr())
  }
  
  # ---- Resolve the shared multivariable settings (defaults for missing) -----
  # Bundled in multi_params so they do not clutter the top-level signature
  mp <- NULL
  if (m$multivariable && !is.null(multi_params)) {
    # NB: screen_tables defaults to NULL (not list()) on purpose — modifyList
    # RECURSES into two lists and merges by name, so an unnamed screen_tables
    # list would be silently dropped against a list() default. NULL forces a
    # straight replacement.
    mp <- utils::modifyList(
      list(
        df = NULL,
        outcome = NULL,
        screen_tables = NULL,
        pval_threshold = 0.2,
        drop_na_cols = NULL
      ),
      multi_params
    )
    # Fail fast (before any fitting) if the screen tables are missing/empty. The
    # count distinguishes an empty list from a list of NULLs (a stale object).
    screen_cols <- colnames(dplyr::bind_rows(mp$screen_tables))
    lacking <- setdiff(c("Predictors", "P-values"), screen_cols)
    if (length(lacking) > 0) {
      stop(
        "compute_all_methods: multi_params$screen_tables lacks ",
        paste(lacking, collapse = ", "),
        " (received ", length(mp$screen_tables), " table(s); columns: ",
        if (length(screen_cols) == 0) "<none>" else paste(screen_cols, collapse = ", "),
        "). Set it to your univariable tables AFTER they are computed.",
        call. = FALSE
      )
    }
  }
  
  # ---- Run the bootstrap INSIDE, once per condition (gated by the toggle) ----
  # Kept in the result so the expensive resampling can be saved to disk
  raw_all <- NULL
  if (m$bootstrap) {
    raw_all <-
      purrr::imap(conditions, function(cond, name) {
        timed(name, "bootstrap", function() {
          bootstrap_lasso_raw(
            cond$x, cond$y,
            B = boot_B, k = boot_k, nfolds = nfolds, alpha = boot_alpha,
            seed = seed
          )
        })
      })
  }
  
  # ---- Compute every enabled block once, per condition (un-prefixed) --------
  # Each condition yields `blocks` (combined-table rows) and `plots` (numeric
  # forest tables), both computed a single time from the same fits
  per_cond <-
    purrr::imap(conditions, function(cond, name) {
      # LASSO selection feeds both lasso_opt AND firth's variable set
      lopt <- if (m$selective_lasso || m$firth) {
        timed(name, "selective LASSO (optimal)", function() {
          si_block(
            cond$x, cond$y,
            k = NULL, alpha = si_alpha, nfolds = nfolds, seed = seed
          )
        })
      } else {
        NULL
      }
      # Accumulate the enabled string blocks (b) and numeric plot tables (p)
      b <- list()
      p <- list()
      # Bootstrap block (summarizes the draws computed above). Guarded: a failed
      # bootstrap leaves NULL draws, and boot_block would error on them.
      if (m$bootstrap && !is.null(raw_all[[name]])) {
        b$boot_lasso <- boot_block(raw_all[[name]], alpha = boot_alpha)
      }
      # Selective inference at the optimal and the min lambda
      if (m$selective_lasso) {
        lmin <- timed(name, "selective LASSO (min)", function() {
          si_block(
            cond$x, cond$y,
            k = 0, alpha = si_alpha, nfolds = nfolds, seed = seed
          )
        })
        b$lasso_opt <- lopt$block
        b$lasso_min <- lmin$block
        p$lasso_opt <- lopt$plot
        p$lasso_min <- lmin$plot
      }
      # Split (multi-sample splitting) block
      if (m$split) {
        sp <- timed(name, "split (cross-fit)", function() {
          split_block(
            cond$x, cond$y,
            n_splits = n_splits, nfolds = nfolds, seed = seed,
            min_pct_chosen = split_min_pct
          )
        })
        b$split <- sp$block
        p$split <- sp$plot
      }
      # Firth block on the LASSO-selected set
      if (m$firth) {
        vars_opt <- if (is.null(lopt)) character() else lopt$block$term
        b$firth <- timed(name, "Firth", function() {
          fit_firth(cond$x, cond$y, vars_opt)
        })
      }
      # Forward stepwise: native (linear selective) + logistic-refit ORs
      if (m$stepwise) {
        sf <- timed(name, "stepwise", function() {
          fit_stepwise(cond$x, cond$y, seed = seed)
        })
        b$step <- sf$logistic
        b$fs <- sf$native
      }
      # Multivariable: computed here per condition (screen univariable, drop this
      # condition's unwanted set, fit one robust logistic). Skipped when the
      # condition sets multi_unwanted = NULL (e.g. collinear or no-surgery).
      if (!is.null(mp) && !is.null(cond$multi_unwanted)) {
        mv <- timed(name, "multivariable", function() {
          fit_multivariable(
            df = mp$df,
            outcome = mp$outcome,
            univ_tables = mp$screen_tables,
            unwanted = cond$multi_unwanted,
            pval_threshold = mp$pval_threshold,
            drop_na_cols = mp$drop_na_cols
          )
        })
        # Guarded: norm_orlike() would error on a NULL table from a failed fit
        if (!is.null(mv)) b$multi <- norm_orlike(mv$table)
      }
      list(blocks = b, plots = p)
    })
  
  # ---- Combined table: prefix each block, join by term, stack ---------------
  # One wide row per (condition, term) with every method's columns side by side
  combined <-
    purrr::imap_dfr(per_cond, function(cc, name) {
      # Empty stub when every method was switched off for this condition
      if (length(cc$blocks) == 0) {
        return(tibble::tibble(term = character(), condition = name))
      }
      # Prefix each block with its method name, then full-join on term
      prefixed <- purrr::imap(cc$blocks, function(tbl, pre) prefix_cols(tbl, pre))
      purrr::reduce(prefixed, dplyr::full_join, by = "term") %>%
        dplyr::mutate(condition = name, .before = 1)
    })
  
  # ---- Attach univariable (the one condition-agnostic method) to "all" -------
  # Univariable tests one predictor at a time, so it is a single full-set table
  if (m$univariable && !is.null(uni_table)) {
    uni_all <-
      norm_orlike(uni_table) %>%
      prefix_cols("uni") %>%
      dplyr::mutate(condition = "all")
    combined <- dplyr::full_join(combined, uni_all, by = c("condition", "term"))
  }
  
  # ---- Sort by selection frequency when the bootstrap ran -------------------
  # Order by bootstrap selection frequency if present, else by term
  if ("boot_lasso_pct_chosen" %in% names(combined)) {
    combined <-
      dplyr::arrange(combined, condition, dplyr::desc(boot_lasso_pct_chosen))
  } else {
    combined <- dplyr::arrange(combined, condition, term)
  }
  
  # ---- Per-method stacked tables (un-prefixed, from the same blocks) --------
  # Stack a given method's block across conditions, tagging the condition
  stack_method <- function(key) {
    purrr::imap_dfr(per_cond, function(cc, name) {
      if (is.null(cc$blocks[[key]])) return(NULL)
      cc$blocks[[key]] %>% dplyr::mutate(condition = name, .before = 1)
    })
  }
  # Assemble the return list: combined first, then each non-empty method table
  results <- list(combined = combined)
  method_keys <- c(
    "lasso_opt", "lasso_min", "boot_lasso", "split", "firth", "step", "fs"
  )
  for (key in method_keys) {
    tbl <- stack_method(key)
    if (nrow(tbl) > 0) results[[key]] <- tbl
  }
  # Univariable is a single full-set table (tagged as the "all" condition)
  if (m$univariable && !is.null(uni_table)) {
    results$uni <-
      norm_orlike(uni_table) %>%
      dplyr::mutate(condition = "all", .before = 1)
  }
  # Keep the bootstrap draws so the caller can saveRDS() them (not a CSV)
  if (!is.null(raw_all)) results$raw_all <- raw_all
  # Keep the numeric forest tables (per condition, per method) for plotting
  results$plot_tables <- purrr::map(per_cond, "plots")
  
  # ---- Report anything that failed ------------------------------------------
  # A failed method drops silently out of the tables, so say so loudly here:
  # a missing column is otherwise indistinguishable from a disabled method.
  results$failures <- failures
  if (verbose && length(failures) > 0) {
    cat(
      "\n", strrep("-", 72), "\n",
      "COMPLETED WITH ", length(failures), " FAILED METHOD(S) — ",
      "their columns are absent from the results:\n",
      paste0("  ", failures, collapse = "\n"), "\n",
      strrep("-", 72), "\n",
      sep = "", file = stderr()
    )
    flush(stderr())
  }
  # Return the named list of tables (nothing written)
  results
}


# ---- The write section (exports only, computes nothing) ---------------------

#' Write a computed-results list to per-method and combined CSVs
#'
#' @description
#' Takes the list returned by [compute_all_methods()] and writes each element to
#' its CSV (rounded for readability). Export only — no computation happens here.
#'
#' @param results Named list from [compute_all_methods()].
#' @param csv_folder Directory to write into.
#' @return `results`, invisibly.
#' @importFrom readr write_csv
#' @importFrom dplyr mutate across where
#' @export
write_all_methods <- function(results, csv_folder) {
  # Map each result element to its output file name
  file_map <- list(
    combined = "results_combined_all_methods.csv",
    lasso_opt = "results_lasso_optimal.csv",
    lasso_min = "results_lasso_min.csv",
    boot_lasso = "results_bootstrap.csv",
    split = "results_split.csv",
    firth = "results_firth.csv",
    step = "results_stepwise_logistic.csv",
    fs = "results_fs_native_linear.csv",
    uni = "results_univariable.csv"
  )
  # OR-scale columns are formatted to display strings so the CSVs match the
  # plots exactly (3 sig figs >= 1 with trailing zeros, 2 dp < 1); every other
  # numeric column is kept numeric at 3 significant figures.
  format_table <- function(df) {
    or_cols <- grep("(^|_)OR($|_)|odds_ratio|Odds Ratios", names(df))
    for (j in or_cols) df[[j]] <- fmt_or_num(df[[j]])
    dplyr::mutate(df, dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 3)))
  }
  # Format and write each element that has a mapped file name
  for (key in names(results)) {
    fname <- file_map[[key]]
    if (is.null(fname)) next
    readr::write_csv(format_table(results[[key]]), file.path(csv_folder, fname))
  }
  # Return the input invisibly for piping
  invisible(results)
}


# ---- The plotting section (forest plots) ------------------------------------

#' Forest plot of odds ratios with CI and p-value text columns
#'
#' @description
#' Renders any table in the shared forest schema (`term`, `odds_ratio`,
#' `pval`, `ci_lo`, `ci_hi`) as a LOG-scale forest plot (point + CI, reference
#' line at OR = 1, point estimate above each dot) alongside right-hand text
#' columns for the confidence interval and the p-value. Used for the
#' selective-inference LASSO (both lambdas) and the split method; the numeric
#' tables it consumes come from `compute_all_methods()$plot_tables`.
#'
#' The axis runs from 0 and its breaks are placed at EQUAL spacing: 0, 0.5, 1,
#' 2, 5, 10 each sit one column apart, with the axis linear in between. It is
#' therefore neither linear (which crushes 0, 0.5 and 1 into the left margin)
#' nor logarithmic (which has no 0). The spacing is a deliberate reading choice,
#' not a property of the numbers: below the null a clinician needs 0, 0.5 and 1
#' and nothing else, and each deserves equal room. 0, 0.5 and 1 are always
#' drawn; the upper ladder is the shortest one that covers the data.
#'
#' @param table A tibble with `term`, `odds_ratio`, `pval`, `ci_lo`, `ci_hi`.
#' @param col_names Named lookup vector passed to [prettify_terms()].
#' @param ci_cap Upper bound of the drawn axis. Interval ends beyond it are
#'   arrowed, and point estimates beyond it are NOT drawn (a clamped marker
#'   would assert a position the estimate does not have) — though the number is
#'   still printed. The CI *text* column always shows the true, uncapped value.
#'   Default 10.
#' @param ci_level Confidence level for the CI-column header (e.g. 0.95 -> "95%
#'   CI"). When `NULL` it is read from the table's `ci_level` attribute (set by
#'   [fit_selective_lasso()] and [split_block()]), falling back to 0.95.
#' @param sub_null_share Fraction of the axis width given to the whole sub-null
#'   block (0 through 0.5 to 1). Raise it to spread the protective end out,
#'   lower it to hand the room to the effects above 1. THIS IS THE ONLY NUMBER
#'   TO TOUCH to retune that balance. Default 0.30.
#' @param pval_digits Decimal places for the p-value column (fixed, so 1 reads
#'   "1.00"); a p below `10^-pval_digits` is shown as e.g. "<0.01".
#' @param widths Relative widths of the (OR plot, spacer, CI, p-value) panels.
#'   The CI panel must stay wide enough for its longest label; the gap between
#'   the two text columns comes from that width alone, since neither column adds
#'   padding of its own.
#'
#' @return A patchwork object (or `NULL`, invisibly, if nothing is estimable),
#'   carrying an `n_terms` attribute so [write_all_plots()] can size the height.
#'
#' @importFrom dplyr filter arrange desc mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_segment geom_vline
#'   geom_text labs guides scale_x_continuous scale_y_discrete expansion theme
#'   element_text element_blank ggtitle theme_void margin coord_cartesian
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom scales trans_new
#' @importFrom grid unit arrow
#' @export
plot_forest <- function(
    table,
    col_names = c(EXPOSURES_CONTINUOUS, EXPOSURES_BINARY, EXPOSURES_CATEGORICAL),
    ci_cap = 10,
    ci_level = NULL,
    sub_null_share = 0.30,
    pval_digits = 2,
    widths = c(1, 0.02, 0.36, 0.14)
) {
  # Nothing to draw for an empty (e.g. nothing-selected) table
  if (is.null(table) || nrow(table) == 0) {
    message("plot_forest: empty table, nothing to plot.")
    return(invisible(NULL))
  }
  # CI level -> header. Read the attribute set upstream (before any dplyr op
  # that would drop it); fall back to the argument, then to 0.95.
  if (is.null(ci_level)) ci_level <- attr(table, "ci_level")
  if (is.null(ci_level)) ci_level <- 0.95
  ci_title <- paste0(round(100 * ci_level), "% CI")
  # Keep only terms with a usable point AND interval (drops split terms chosen
  # in 0-1 splits: a forest entry needs both to be meaningful)
  coef_table <-
    table %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::filter(is.finite(odds_ratio), is.finite(ci_lo), is.finite(ci_hi))
  # Nothing left to draw once un-estimable terms are removed
  if (nrow(coef_table) == 0) {
    message("plot_forest: no estimable terms, nothing to plot.")
    return(invisible(NULL))
  }
  # Row order: with more than 3 predictors, sort by p-value (smallest at the
  # top); with 3 or fewer, sort by effect size. ggplot draws the first factor
  # level at the bottom, so the top row must be the LAST level -> descending.
  coef_table <- if (nrow(coef_table) > 3) {
    dplyr::arrange(coef_table, dplyr::desc(pval))
  } else {
    dplyr::arrange(coef_table, dplyr::desc(odds_ratio))
  }
  coef_table <-
    coef_table %>%
    dplyr::mutate(
      term = prettify_terms(term, col_names),
      term = factor(term, levels = term)
    ) %>%
    dplyr::mutate(
      # Keep the true values: the text columns are never capped
      or_raw = odds_ratio,
      ci_hi_raw = ci_hi,
      # Flag whatever runs past the top of the axis BEFORE clamping
      or_off = odds_ratio > ci_cap,
      hi_off = ci_hi > ci_cap,
      # Clamp to the drawn range; the axis is bounded below by 0 anyway
      odds_ratio = pmin(odds_ratio, ci_cap),
      ci_hi = pmin(ci_hi, ci_cap)
    )
  # OR / CI number labels use the shared OR formatter, so the plot matches the
  # tables exactly (3 sig figs >= 1 with trailing zeros, 2 dp < 1)
  fmt_num <- fmt_or_num
  # P-value formatting, by magnitude:
  #   >= 1e-2  : 2 significant figures, min 2 dp (so p = 1 reads "1.00")
  #   1e-4..1e-2: 1 significant figure (0.001 -> "0.001", 0.00068 -> "0.0007")
  #   < 1e-4   : scientific, 1 significant figure (e.g. "7e-05")
  fmt_p <- function(p) {
    vapply(p, function(x) {
      if (is.na(x)) return(NA_character_)
      if (x < 1e-4) return(formatC(x, format = "e", digits = 0))
      if (x < 1e-2) {
        # Decimals from the ROUNDED value so 0.001 is "0.001", not "0.0010"
        r <- signif(x, 1)
        return(formatC(r, format = "f", digits = -floor(log10(r))))
      }
      formatC(x, format = "f", digits = max(2, 1 - floor(log10(x))))
    }, character(1))
  }
  # Breaks. 0, 0.5 and 1 are ALWAYS drawn: below the null those are the only
  # values a reader needs, and 1 must be visible to locate the null. The upper
  # ladder is the shortest one that still covers the data.
  data_max <- min(max(coef_table$ci_hi, coef_table$odds_ratio), ci_cap)
  ladders <- list(
    c(0, 0.5, 1, 2),
    c(0, 0.5, 1, 2, 3, 4),
    c(0, 0.5, 1, 2, 4, 6),
    c(0, 0.5, 1, 2, 4, 6, 8),
    c(0, 0.5, 1, 2, 5, 10)
  )
  idx <- which(vapply(ladders, function(b) max(b) >= data_max, logical(1)))[1]
  x_breaks <- if (is.na(idx)) c(0, 0.5, 1, 2, 5, ci_cap) else ladders[[idx]]
  x_hi <- max(x_breaks)
  # Trailing zeros are noise on an axis: "0.5", not "0.50"; "10", not "10.0".
  # trim = TRUE is essential: format() otherwise pads every label to a COMMON
  # width, and the padding shifts each label off the center of its gridline.
  x_labels <- format(x_breaks, trim = TRUE, drop0trailing = TRUE)
  # Display transform: the break positions are assigned by us, not derived from
  # the numbers, with the axis linear between them. No analytic curve can do
  # this — a ladder unevenly spaced in value cannot come out evenly spaced under
  # linear, log, or power.
  # ===========================================================================
  # >>> THE KNOB FOR HOW MUCH WIDTH THE 0-0.5-1 BLOCK GETS IS `sub_null_share`,
  # >>> IN THIS FUNCTION'S ARGUMENTS ABOVE. Nothing here needs editing.
  # ===========================================================================
  # Slot widths. Every interval at or below the null gets one unit; every
  # interval above it gets `upper_weight`, solved so that the sub-null block
  # ends up occupying exactly `sub_null_share` of the axis. Within each group
  # the slots stay equal, so 0-0.5 = 0.5-1 and 1-2 = 2-5 = 5-10 as required.
  n_sub <- sum(x_breaks[-1] <= 1)
  n_upper <- sum(x_breaks[-1] > 1)
  upper_weight <- (n_sub * (1 - sub_null_share)) / (n_upper * sub_null_share)
  knots <- x_breaks
  slots <- c(0, cumsum(ifelse(x_breaks[-1] <= 1, 1, upper_weight)))
  x_trans <- scales::trans_new(
    name = "or_slots",
    transform = function(x) stats::approx(knots, slots, xout = x, rule = 2)$y,
    inverse = function(x) stats::approx(slots, knots, xout = x, rule = 2)$y
  )
  # Tight expansion so gridlines hug the data; clip = "off" lets the OR labels
  # sit above the top dot without being cut. Outer space comes from plot.margin.
  y_expand <- ggplot2::expansion(add = c(0.15, 0.2))
  pm <- grid::unit(c(1.4, 0, 1.4, 0), "lines")
  # CI column: a left margin sets it off from the plot
  pm_ci <- grid::unit(c(1.4, 0, 1.4, 0.8), "lines")
  # P-value column: NO left margin. Its separation from the CI column comes from
  # the CI column's own width, so no dead space is added on top of that. The
  # right margin just keeps the last column off the edge of the figure.
  pm_pval <- grid::unit(c(1.4, 1, 1.4, 0), "lines")
  # Arrowheads sit a fixed step inside the axis end
  arrow_head <- grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")
  arrow_step <- 0.97
  odds_ratios_plot <-
    coef_table %>%
    ggplot2::ggplot(ggplot2::aes(x = odds_ratio, y = term, color = term)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
    # The interval, clamped to the axis
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_lo, xend = ci_hi, y = term, yend = term)
    ) +
    # A whisker tick marks the lower end, and the upper end when it fits
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lo, xmax = ci_lo), width = 0.08
    ) +
    ggplot2::geom_errorbar(
      data = function(d) dplyr::filter(d, !hi_off),
      ggplot2::aes(xmin = ci_hi, xmax = ci_hi), width = 0.08
    ) +
    # An upper end that runs past the axis gets an arrowhead instead of a tick,
    # so it reads as "continues beyond" rather than "ends here"
    ggplot2::geom_segment(
      data = function(d) dplyr::filter(d, hi_off),
      ggplot2::aes(x = ci_hi * arrow_step, xend = ci_hi, y = term, yend = term),
      arrow = arrow_head
    ) +
    # The marker is drawn ONLY when the estimate is on the axis. A clamped dot
    # would place the point somewhere it is not; better to show no dot at all.
    ggplot2::geom_point(data = function(d) dplyr::filter(d, !or_off)) +
    # The number, however, is always printed (at the axis edge when off-scale),
    # so an unplottable estimate is still reported rather than hidden
    ggplot2::geom_text(
      ggplot2::aes(label = fmt_num(or_raw)),
      vjust = -0.9, size = 3.3, show.legend = FALSE
    ) +
    ggplot2::labs(x = "Odds Ratio", y = NULL) +
    ggplot2::guides(color = "none") +
    ggplot2::scale_x_continuous(
      trans = x_trans, breaks = x_breaks, labels = x_labels,
      limits = c(0, x_hi),
      expand = ggplot2::expansion(mult = c(0.02, 0.03))
    ) +
    ggplot2::scale_y_discrete(limits = levels(coef_table$term), expand = y_expand) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      # hjust = 0.5 pins each label to the center of its own gridline
      axis.text.x = ggplot2::element_text(
        size = 12, hjust = 0.5, margin = ggplot2::margin(t = 6)
      ),
      axis.text.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 6)),
      axis.title = ggplot2::element_text(size = 13),
      plot.margin = pm
    )
  # Values start hard at the panel's left edge (x = 0, hjust = 0) with NO
  # expansion on either side, so the column adds no dead space of its own: its
  # width alone decides the gap to the column that follows.
  x_left <- ggplot2::scale_x_continuous(
    limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0))
  )
  # Header is LEFT-aligned to the same panel-left edge as its values (hjust = 0),
  # so the two share a left edge no matter how wide the column is. Margin lifts
  # it off the first row.
  title_theme <- ggplot2::element_text(
    hjust = 0, face = "bold", margin = ggplot2::margin(b = 10)
  )
  ci_plot <-
    coef_table %>%
    dplyr::mutate(
      # The RAW upper bound: the text column reports the interval, uncapped
      ci_label = paste0(fmt_num(ci_lo), " – ", fmt_num(ci_hi_raw))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = 0, y = term)) +
    ggplot2::geom_text(ggplot2::aes(label = ci_label), hjust = 0, size = 4) +
    x_left +
    ggplot2::scale_y_discrete(limits = levels(coef_table$term), expand = y_expand) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(ci_title) +
    ggplot2::theme(plot.margin = pm_ci, plot.title = title_theme)
  pvalues_table_plot <-
    coef_table %>%
    dplyr::mutate(pval = fmt_p(pval)) %>%
    ggplot2::ggplot(ggplot2::aes(x = 0, y = term)) +
    ggplot2::geom_text(ggplot2::aes(label = pval), hjust = 0, size = 4) +
    x_left +
    ggplot2::scale_y_discrete(limits = levels(coef_table$term), expand = y_expand) +
    ggplot2::theme_void() +
    ggplot2::ggtitle("P-value") +
    ggplot2::theme(plot.margin = pm_pval, plot.title = title_theme)
  combined_plot <-
    odds_ratios_plot +
    patchwork::plot_spacer() +
    ci_plot +
    pvalues_table_plot +
    patchwork::plot_layout(ncol = 4, widths = widths)
  # Carry the predictor count so write_all_plots can size the figure height
  attr(combined_plot, "n_terms") <- nrow(coef_table)
  combined_plot
}


#' Build a forest plot for every selective-LASSO (both lambdas) and split table
#'
#' @description
#' Maps over `compute_all_methods()$plot_tables` and builds one [plot_forest()]
#' per (condition, method) for the methods that carry a numeric OR table:
#' `lasso_opt`, `lasso_min` and `split`. Builds objects only — writing is left
#' to [write_all_plots()]. Empty (nothing-selected) tables are skipped.
#'
#' @param results Named list from [compute_all_methods()] (uses `plot_tables`).
#' @param col_names Named lookup vector passed through to [plot_forest()].
#' @param methods Named logical list toggling which plottable methods to render
#'   (`lasso_opt`, `lasso_min`, `split`); missing entries default to `TRUE`.
#' @param ... Further arguments forwarded to [plot_forest()] (e.g. `ci_cap`).
#' @return Named list of patchwork objects, keyed `"<label>_<condition>"`.
#' @importFrom purrr map
#' @export
plot_all_methods <- function(
    results,
    col_names = c(EXPOSURES_CONTINUOUS, EXPOSURES_BINARY, EXPOSURES_CATEGORICAL),
    methods = NULL,
    ...
) {
  # Human-readable labels used in the returned names (and hence file names)
  plot_labels <- c(
    lasso_opt = "lasso_optimal",
    lasso_min = "lasso_min",
    split = "split"
  )
  # Resolve which methods to plot (any missing entry defaults to TRUE)
  m <- utils::modifyList(
    list(lasso_opt = TRUE, lasso_min = TRUE, split = TRUE),
    if (is.null(methods)) list() else methods
  )
  # Nothing to do if compute_all_methods kept no forest tables
  pt <- results$plot_tables
  if (is.null(pt)) return(list())
  # One forest per (condition, method), skipping empty/unselected tables
  out <- list()
  for (cond in names(pt)) {
    for (method in names(pt[[cond]])) {
      # Skip methods switched off for plotting
      if (!isTRUE(m[[method]])) next
      tbl <- pt[[cond]][[method]]
      if (is.null(tbl) || nrow(tbl) == 0) next
      # Fall back to the raw method name if it has no pretty label
      label <- plot_labels[[method]]
      if (is.null(label)) label <- method
      key <- paste(label, cond, sep = "_")
      out[[key]] <- plot_forest(tbl, col_names = col_names, ...)
    }
  }
  out
}


#' Write forest plots to file (exports only, builds nothing)
#'
#' @description
#' Takes the list from [plot_all_methods()] and saves each as
#' `forest_<label>_<condition>.<format>`. Export only — no plotting logic here.
#' Defaults to PDF: forest plots are line art, so vector output is what a journal
#' wants, and it stays sharp at any size.
#'
#' @param plots Named list from [plot_all_methods()].
#' @param plot_folder Directory to write into.
#' @param formats Character vector of output formats, e.g. `"pdf"` (default),
#'   `"png"`, or `c("pdf", "png")` to write both.
#' @param width Figure width in inches.
#' @param height Figure height in inches; `NULL` (default) auto-sizes by the
#'   predictor count in three tiers (3 / 6 / 10 in for 1-3 / 4-10 / >10 terms)
#'   so sparse plots are not mostly empty space. Pass a number to fix it.
#' @param dpi Resolution (raster formats only; ignored for PDF).
#' @return `plots`, invisibly.
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices cairo_pdf
#' @export
write_all_plots <- function(
    plots,
    plot_folder,
    formats = "pdf",
    width = 9,
    height = NULL,
    dpi = 300
) {
  # Three-tier figure height by predictor count (used when height is not fixed)
  tier_height <- function(n) {
    if (is.null(n)) return(6)
    if (n <= 3) 3 else if (n <= 10) 6 else 10
  }
  # cairo_pdf, not the base pdf device: the CI labels contain a Unicode en-dash,
  # which the base device silently drops
  device_for <- function(fmt) {
    if (identical(fmt, "pdf")) grDevices::cairo_pdf else fmt
  }
  # Save every non-NULL forest under its already-file-friendly key
  for (key in names(plots)) {
    p <- plots[[key]]
    if (is.null(p)) next
    # Auto-size height from the predictor count unless the caller fixed it
    h <- if (!is.null(height)) height else tier_height(attr(p, "n_terms"))
    for (fmt in formats) {
      ggplot2::ggsave(
        filename = file.path(plot_folder, paste0("forest_", key, ".", fmt)),
        plot = p,
        device = device_for(fmt),
        width = width,
        height = h,
        dpi = dpi
      )
    }
  }
  # Return the input invisibly for piping
  invisible(plots)
}