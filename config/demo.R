# Settings for one run. Copy this file for a new export and run
#   Rscript R/run_all.R config/<your file>.R
# Everything a lab needs to change is here; the scripts in R/ should not need editing.

cfg <- list(
  name = "Synthetic test data (AMR package example_isolates)",
  demo = TRUE,                      # TRUE only for the bundled example: enables the hand-written
                                    # interpretation in R/03_report.R and the recovery test
  input = "data/raw/isolates_messy.csv",
  column_map = NULL,                # CSV mapping export columns to standard names; NULL if the
                                    # export already uses them (see config/column_map_whonet.csv)
  salt_file = "salt.txt",           # secret used to pseudonymise patient IDs; never commit it
  processed = "data/processed/isolates_clean.csv",
  reports = "reports/demo",

  # Breakpoints for raw MIC and disk values. S/I/R text from the lab is kept as reported.
  guideline = "EUCAST 2025",

  first_isolate_method = "episode-based",
  episode_days = 365,               # first isolate per patient and species per 365 days (CLSI M39)
  years = 2002:2017,                # analysis window; NULL = every year in the data
  antibiogram_years = 2015:2017,    # the printed antibiogram; CLSI M39 asks for the latest year,
                                    # pooled with earlier ones only when n < 30
  period_length = 4,                # years per period in the trend table and chart
  minimum = 30,                     # no percentage is reported for fewer tested isolates
  # Gram-negative and Gram-positive panels are kept apart so that drugs an organism is
  # intrinsically resistant to (e.g. vancomycin for E. coli) do not show up as "resistance".
  antibiogram_gram_negative = c("AMX", "AMC", "CXM", "CTX", "CAZ", "MEM", "GEN", "CIP", "SXT", "NIT"),
  antibiogram_gram_positive = c("PEN", "FLC", "ERY", "CLI", "GEN", "CIP", "SXT", "TCY", "VAN", "RIF"),
  pairs = data.frame(mo = "B_ESCHR_COLI", ab = c("CIP", "AMC", "SXT"))
)
