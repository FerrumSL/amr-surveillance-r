# Export the AMR package's bundled example data to data/raw/ and add deliberate mess.
#
# The source is `AMR::example_isolates`: 2,000 randomised, fictitious isolates that the
# package authors describe as reflecting reality. It is NOT real surveillance data.
# The mess below is introduced on purpose, with a fixed seed, so that R/01_clean.R has
# realistic problems to fix (and so the test can check that cleaning recovers the original).

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
})

set.seed(2026)

# Replace a random share of values in x with `variants(x)`.
messify <- function(x, share, variants) {
  hit <- !is.na(x) & runif(length(x)) < share
  x[hit] <- variants(x[hit])
  x
}

sir_variants <- function(v) {
  out <- v
  for (code in c("S", "I", "R")) {
    idx <- which(v == code)
    choices <- switch(code,
      S = c("s", " S", "S ", "Susceptible", "SENSITIVE"),
      I = c("i", "Intermediate", "I "),
      R = c("r", "R ", "Resistant", "RESISTANT")
    )
    out[idx] <- sample(choices, length(idx), replace = TRUE)
  }
  out
}

org_variants <- function(mo_code) {
  name <- mo_name(mo_code, language = NULL)
  style <- sample(1:4, length(mo_code), replace = TRUE)
  ifelse(style == 1, toupper(name),
    ifelse(style == 2, paste0(" ", tolower(name), " "),
      ifelse(style == 3, mo_shortname(mo_code, language = NULL), name)))
}

ab_cols <- names(example_isolates)[vapply(example_isolates, is.sir, logical(1))]

raw <- example_isolates |>
  mutate(
    organism = as.character(mo),
    across(all_of(ab_cols), as.character),
    date = format(date, "%Y-%m-%d"),
    gender = as.character(gender)
  ) |>
  select(-mo) |>
  relocate(organism, .after = ward)

raw$organism <- messify(raw$organism, 0.6, org_variants)
raw[ab_cols] <- lapply(raw[ab_cols], messify, share = 0.15, variants = sir_variants)
# Some labs export day-first dates; this is the only other format the cleaner accepts.
raw$date <- messify(raw$date, 0.1, function(d) format(as.Date(d), "%d/%m/%Y"))
raw$gender <- messify(raw$gender, 0.2, tolower)
raw$ward <- messify(raw$ward, 0.1, function(w) paste0(w, " "))

write_csv(raw, "data/raw/isolates_messy.csv", na = "")
message("Wrote data/raw/isolates_messy.csv (", nrow(raw), " rows)")
