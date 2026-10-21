# Run the whole workflow from the project root:
#   Rscript R/run_all.R                    demo data, settings in config/demo.R
#   Rscript R/run_all.R config/ncbi.R      real NCBI isolates
#   Rscript R/run_all.R config/my_lab.R    a lab's own export (copy config/demo.R)
# The demo and NCBI inputs are downloaded or generated only if missing, so a lab's own
# export can be dropped in their place.

args <- commandArgs(trailingOnly = TRUE)
source(if (length(args) > 0) args[1] else "config/demo.R")

if (cfg$input == "data/raw/isolates_messy.csv" && !file.exists(cfg$input)) source("R/00_fetch.R")
if (cfg$input == "data/raw/ncbi_export.csv" && !file.exists(cfg$input)) source("R/00_fetch_ncbi.R")
source("R/01_clean.R")
source("R/02_analyse.R")
source("R/03_report.R")
if (cfg$input == "data/raw/ncbi_export.csv") source("R/04_validate.R")
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
