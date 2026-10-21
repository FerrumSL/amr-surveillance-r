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
  expect_equal(n_distinct(got$patient), n_distinct(want$patient))
  for (ab in names(want)[vapply(want, is.sir, logical(1))]) {
    expect_equal(as.character(got[[ab]]), as.character(want[[ab]]), label = ab)
  }
})

test_that("antibiotic columns are matched exactly, never guessed", {
  got <- antibiotic_columns(c("CIP", "Ciprofloxacin", "CIP_ND5", "AMC_EM", "PEN_NE", "SXT_ND1_2",
                              "Country", "Ward", "ESBL", "comments"))
  expect_equal(got$ab, c("CIP", "CIP", "CIP", "AMC", "PEN", "SXT", NA, NA, NA, NA))
  expect_equal(got$method, c(NA, NA, "disk", "mic", "mic", "disk", NA, NA, NA, NA))
  expect_equal(got$system, c(NA, NA, "CLSI", "EUCAST", "CLSI", "CLSI", NA, NA, NA, NA))
})

test_that("raw MIC and disk values are interpreted; S/I/R text is kept as reported", {
  mo <- as.mo(rep("E. coli", 5), info = FALSE)
  mic <- interpret_column(c("<=0.25", ">4", "R", "", "16/4"), mo, "CIP", "mic", "EUCAST 2025")
  expect_equal(as.character(mic$sir), c("S", "R", "R", NA, "R"))
  expect_equal(mic$n_numeric, 3)
  disk <- interpret_column(c("30", "10"), mo[1:2], "CIP", "disk", "EUCAST 2025")
  expect_equal(as.character(disk$sir), c("S", "R"))
  expect_error(interpret_column(c("0.5"), mo[1], "CIP", NA, "EUCAST 2025", "Cipro"), "Cipro")
  expect_error(interpret_column(c("0,5"), mo[1], "CIP", NA, "EUCAST 2025", "Cipro"), "Cipro")
  # Lab spellings: decimal comma, "=<", units. Low (susceptible) MICs must not be lost.
  spelled <- interpret_column(c("0,25", "=<0.25", "0.25 mg/L", ">=6.4/121.6x"), mo[1:4], "CIP", "mic", "EUCAST 2025")
  expect_equal(as.character(spelled$sir), c("S", "S", "S", NA))
  expect_equal(spelled$invalid, ">=6.4/121.6x")
  expect_equal(as.character(interpret_column("30 mm", mo[1], "CIP", "disk", "EUCAST 2025")$sir), "S")
})

test_that("patient IDs are replaced by stable keyed hashes", {
  salt <- strrep("a", 32)
  h <- pseudonymise(c("MRN001", "MRN002", "MRN001", NA), salt)
  expect_equal(h[1], h[3])
  expect_false(h[1] == h[2])
  expect_true(is.na(h[4]))
  expect_false(h[1] == pseudonymise("MRN001", strrep("b", 32)))
  raw <- tibble(date = "2015-10-06", patient = "MRN001", organism = "E. coli", CIP = "S")
  expect_false("MRN001" %in% unlist(clean_isolates(raw, lookup)$data))
})

test_that("a WHONET export loads through the WHONET column map", {
  raw <- mutate(WHONET, across(everything(), as.character))
  map <- read_csv(file.path(root, "config/column_map_whonet.csv"), col_types = cols(.default = col_character()))
  res <- clean_isolates(raw, lookup, map)
  d <- res$data
  expect_equal(nrow(d), nrow(distinct(d)))
  expect_false(any(c("Last name", "First name") %in% names(d)))
  expect_false(any(raw$`Last name` %in% unlist(d)))
  expect_true(all(c("AMP", "CIP", "CLI") %in% names(d)))  # CIP_ED5 and CIP_EE merged into CIP
  icr <- toupper(raw$`Inducible clindamycin resistance`) == "TRUE"
  expect_true(all(d$CLI[icr[seq_len(nrow(d))]] == "R", na.rm = TRUE))
})

test_that("alerts catch MRSA, VRE and CRE but not intrinsic resistance", {
  d <- tibble(date = as.Date("2020-01-01"), patient = c("a", "b", "c", "d", "e"), ward = "ICU",
              specimen = NA_character_,
              mo = as.mo(c("S. aureus", "E. faecium", "E. casseliflavus", "K. pneumoniae", "P. mirabilis"),
                         info = FALSE),
              FOX = as.sir(c("R", NA, NA, NA, NA)), VAN = as.sir(c(NA, "R", "R", NA, NA)),
              IPM = as.sir(c(NA, NA, NA, "S", "R")), MEM = as.sir(c(NA, NA, NA, "R", "S")))
  got <- flag_alerts(d)
  expect_equal(sort(paste(got$alert, got$patient)), c("CRE d", "MRSA a", "VRE b"))
  # Unmatched organisms never alert, and no alerts gives an empty table, not an error.
  none <- mutate(d[1, ], mo = as.mo(NA, info = FALSE), FOX = as.sir("S"), VAN = as.sir("R"))
  expect_equal(nrow(flag_alerts(none)), 0)
  expect_true("alert" %in% names(flag_alerts(none)))
})

test_that("zones from another system's disks are not interpreted; MICs always are", {
  raw <- tibble(date = "2020-01-01", patient = "1", organism = "E. coli",
                CIP_ND5 = "20", CIP_ED5 = "30", AMC_NM = "4", AMP_BD10 = "25")
  res <- clean_isolates(raw, lookup, guideline = "EUCAST 2025")
  expect_equal(res$log$antibiotic_columns$guideline, c(NA, "EUCAST 2025", "EUCAST 2025", NA))
  expect_equal(as.character(res$data$CIP), "S")   # from CIP_ED5, not the CLSI zone
  expect_true(is.na(res$data$AMP))
  expect_equal(sum(res$log$skipped$n), 2)
})
