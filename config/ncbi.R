# Real AST results from NCBI BioSample (human isolates), made by R/00_fetch_ncbi.R.
# Run with: Rscript R/run_all.R config/ncbi.R
#
# Not surveillance data. Submitters chose what to sequence, so resistance is over-represented,
# and each isolate is treated as a separate patient because NCBI has no patient IDs.

cfg <- list(
  name = "Example: public NCBI BioSample isolates (not surveillance data)",
  demo = FALSE,
  input = "data/raw/ncbi_export.csv",
  column_map = "config/column_map_ncbi.csv",
  salt_file = "salt.txt",
  processed = "data/processed/ncbi_clean.csv",
  reports = "reports/ncbi",

  # Most submitters used CLSI. Every raw value is re-interpreted with one version, as
  # surveillance should be; R/04_validate.R measures how often that differs from the lab call.
  guideline = "CLSI 2025",

  first_isolate_method = "episode-based",
  episode_days = 365,
  years = NULL,
  antibiogram_years = NULL,         # NULL = all years (no single year has enough isolates)
  period_length = 5,
  minimum = 30,
  antibiogram_gram_negative = c("AMP", "AMC", "TZP", "CRO", "CAZ", "FEP", "MEM", "IPM", "GEN",
                                "AMK", "CIP", "SXT", "COL"),
  antibiogram_gram_positive = c("PEN", "OXA", "ERY", "CLI", "GEN", "CIP", "SXT", "TCY", "VAN",
                                "LNZ", "DAP", "RIF"),
  pairs = data.frame(mo = c("B_KLBSL_PNMN", "B_ESCHR_COLI"), ab = c("MEM", "CIP")),

  # Printed on the antibiogram (optional in any config).
  antibiogram_subtitle = "All years pooled: no single year has 30 first isolates per organism",
  antibiogram_note = paste("Public NCBI submissions, enriched for resistant isolates.",
                           "Not representative of any population and not for empirical therapy decisions.")
)
