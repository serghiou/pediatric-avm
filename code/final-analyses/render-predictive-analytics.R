# =============================================================================
# Render the predictive-analytics notebook once per outcome.
# -----------------------------------------------------------------------------
# The heavy methods (bootstrap, cross-fitted split) parallelize through furrr.
# The future plan is set ONCE here and reused by every render: setting it inside
# the .Rmd would spin up and tear down a fresh set of worker sessions for every
# outcome. The plan lives in the future package's own state, not in an
# environment, so each render still gets a clean `new.env()` of its own.
# =============================================================================

library(rmarkdown)
library(furrr)


# ---- Parameters -------------------------------------------------------------

proportion_of_data <- 1

rmd <- "code/final-analyses/predictive-analytics.Rmd"

configs <- c(
  "config_is-poor-outcome.R",
  "config_is-poor-outcome_haemTRUE.R",
  "config_is-poor-outcome_haemFALSE.R",
  "config_has-complication-major.R",
  "config_has-complication-minor.R",
  "config_is-obliterated.R"
)

titles <- c(
  "Poor outcome",
  "Poor outcome - Hemorrhage",
  "Poor outcome - No Hemorrhage",
  "Complications - Major",
  "Complications - Minor",
  "Obliterated"
)

# Fail loudly now rather than mislabel an outcome later
stopifnot(length(configs) == length(titles))


# ---- Parallel backend -------------------------------------------------------

# Leave one core for the OS; multisession works on macOS, Windows and Linux
future::plan(future::multisession, workers = future::availableCores() - 1)


# ---- Render -----------------------------------------------------------------

# Collect failures so one bad config does not abort the whole batch
failures <- character()

for (i in seq_along(configs)) {
  
  params <- list(
    PROPORTION_OF_DATA = proportion_of_data,
    CONFIG = configs[i],
    TITLE = titles[i]
  )
  
  message("→ Rendering: ", titles[i])
  started <- Sys.time()
  
  # A fresh environment per outcome. Without it, a config that fails to define
  # something would silently inherit the PREVIOUS outcome's value — a stale
  # UNI_TABLE or unwanted_* set would corrupt the results without a word.
  env <- new.env()
  
  tryCatch({
    out <- rmarkdown::render(rmd, params = params, envir = env)
    
    # The notebook defines `html_folder` (dated) at its end; file the HTML there
    html_dir <- file.path(dirname(rmd), get("html_folder", envir = env))
    dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Filename = config name without the "config_" prefix or the ".R" extension
    stub <- sub("^config_", "", tools::file_path_sans_ext(configs[i]))
    dest <- file.path(html_dir, paste0(stub, ".html"))
    file.copy(out, dest, overwrite = TRUE)
    file.remove(out)
    
    mins <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)
    message("✔ Done: ", titles[i], "  →  ", dest, "  (", mins, " min)")
  },
  error = function(e) {
    message("✖ FAILED: ", titles[i], " - ", conditionMessage(e))
    failures <<- c(failures, configs[i])
  })
}


# ---- Report -----------------------------------------------------------------

if (length(failures) == 0) {
  message("All ", length(configs), " analyses complete.")
} else {
  message(
    "Completed ", length(configs) - length(failures), " / ", length(configs),
    ". Failed: ", paste(failures, collapse = ", ")
  )
}


# ---- Hand the cores back ----------------------------------------------------

future::plan(future::sequential)