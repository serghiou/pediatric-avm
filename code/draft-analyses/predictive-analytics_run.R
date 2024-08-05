# Load libraries
library(rmarkdown)
library(stringr)

# Define parameters
proportion_of_data <- 1

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

# Create the output file names
output_files_prefix <- "../../outputs/predictive-analytics/"
output_files_suffix <- stringr::str_extract(configs, "(?<=config_)[^.]+(?=\\.R)")
output_files <- paste0(output_files_prefix, output_files_suffix, ".html")

# Render
for (i in seq_along(configs)) {
  # Create parameters list
  params <- list(
    PROPORTION_OF_DATA = proportion_of_data,
    CONFIG = configs[i],
    TITLE = titles[i]
  )
  
  # Render
  rmarkdown::render(
    input = "code/draft-analyses/predictive-analytics.Rmd",
    params = params,
    output_file = output_files[i],
    envir = new.env()
  )
}
