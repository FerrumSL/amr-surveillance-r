# Run from the project root: Rscript -e 'testthat::test_dir("tests/testthat")'
suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
})
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
withr::local_options(amrwf.skip_main = TRUE)
source(file.path(root, "R/01_clean.R"), local = TRUE)
lookup <- read_csv(file.path(root, "data/organism_lookup.csv"),
                   col_types = cols(.default = col_character()))

test_that("messy S/I/R text maps to the right category", {
  x <- c("S", "s", " S", "Susceptible", "SENSITIVE", "i", "Intermediate", "r", "R ",
         "Resistant", NA, "")
  expect_equal(as.character(clean_sir(x)),
               c(rep("S", 5), "I", "I", "R", "R", "R", NA, NA))
})

test_that("ISO and dd/mm/yyyy dates parse; anything else is NA", {
  expect_equal(clean_date(c("2015-10-06", "06/10/2015", " 2015-10-06 ", "10/6/15", "10-06-2015",
                            "15-10-06", "2015-10-06 junk", "31/02/2015")),
               as.Date(c("2015-10-06", "2015-10-06", "2015-10-06", NA, NA, NA, NA, NA)))
})

test_that("columns that are not exact antibiotic codes stop the run", {
  raw <- tibble(date = "2015-10-06", patient = "1", age = "50", gender = "F", ward = "ICU",
                organism = "E. coli", CIP = "S", comments = "S")
  expect_error(clean_isolates(raw, lookup), "comments")
})

test_that("ambiguous abbreviations follow the reviewed lookup", {
  got <- as.character(clean_organism(c("S. hominis", " s. HOMINIS ", "P. vulneris", "E. coli",
                                       "STAPHYLOCOCCUS AUREUS"), lookup)$mo)
  expect_equal(got, c("B_STPHY_HMNS", "B_STPHY_HMNS", "B_PSDSC_VLNR", "B_ESCHR_COLI",
                      "B_STPHY_AURS"))
})

test_that("every organism in the example data survives every spelling used in data/raw", {
  mos <- unique(as.character(example_isolates$mo))
  for (spelling in list(mo_name(mos, language = NULL), toupper(mo_name(mos, language = NULL)),
                        tolower(mo_name(mos, language = NULL)),
                        mo_shortname(mos, language = NULL), mos)) {
    expect_equal(as.character(clean_organism(spelling, lookup)$mo), mos)
  }
})

test_that("cleaning the messy export recovers the package's original data exactly", {
  raw <- read_csv(file.path(root, "data/raw/isolates_messy.csv"),
                  col_types = cols(.default = col_character()), na = "")
  skip_if_not(identical(trimws(raw$patient), as.character(example_isolates$patient)),
              "data/raw/isolates_messy.csv is not the demo export")
  got <- clean_isolates(raw, lookup)$data
  want <- distinct(example_isolates)
  expect_equal(nrow(got), nrow(want))
  expect_equal(got$date, want$date)
  expect_equal(as.character(got$mo), as.character(want$mo))
  for (ab in names(want)[vapply(want, is.sir, logical(1))]) {
    expect_equal(as.character(got[[ab]]), as.character(want[[ab]]), label = ab)
  }
})
