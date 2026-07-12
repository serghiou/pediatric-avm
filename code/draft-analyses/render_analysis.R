# Renders the notebook for a given config + title, then files the produced
# HTML under the html/ folder that the notebook itself creates at the end.
render_one <- function(
    config,
    title,
    proportion = 1,
    rmd = "code/draft-analyses/predictive-analytics.Rmd") {
  
  message("→ Rendering: ", config, "  (", title, ")")
  
  out <- rmarkdown::render(
    input = rmd,
    output_format = "rmdformats::readthedown",
    params = list(
      CONFIG = config,
      TITLE = title,
      PROPORTION_OF_DATA = proportion
    ),
    envir = globalenv()
  )
  
  # The notebook defines `html_folder` near the end; reuse it
  html_folder <- get("html_folder", envir = globalenv())
  html_dir <- normalizePath(file.path(dirname(rmd), html_folder), mustWork = FALSE)
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Filename = config name without "config_" prefix or ".R"
  stub <- sub("^config_", "", tools::file_path_sans_ext(config))
  dest <- file.path(html_dir, paste0(stub, ".html"))
  file.copy(out, dest, overwrite = TRUE)
  file.remove(out)  # remove the copy left next to the .Rmd
  
  message("✔ Done: ", config, "  →  ", dest)
  invisible(dest)
}

# --- Run analysis ---
# To only run one config:
# render_one("config_has-complication-minor.R", "Complications - Minor")

runs <- tibble::tribble(
  ~config,                             ~title,
  "config_is-poor-outcome.R",          "Poor outcome",
  "config_is-poor-outcome_haemTRUE.R", "Poor outcome - Hemorrhage",
  "config_is-poor-outcome_haemFALSE.R","Poor outcome - No Hemorrhage",
  "config_has-complication-major.R",   "Complications - Major",
  "config_has-complication-minor.R",   "Complications - Minor",
  "config_is-obliterated.R",           "Obliterated"
)

purrr::pwalk(runs, \(config, title) render_one(config, title))
message("All ", nrow(runs), " analyses complete.")
