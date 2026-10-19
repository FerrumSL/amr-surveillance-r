# Run the whole workflow from the project root: Rscript R/run_all.R
# data/raw/isolates_messy.csv is created by R/00_fetch.R only if it is missing, so a lab's
# own export can be dropped in its place.

if (!file.exists("data/raw/isolates_messy.csv")) source("R/00_fetch.R")
source("R/01_clean.R")
source("R/02_analyse.R")
source("R/03_report.R")
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
