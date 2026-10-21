# Clean, structure and pseudonymise a raw isolate export, and list infection-control alerts.
#
# Input:  cfg$input       one row per isolate; any column names, mapped by cfg$column_map
#         cfg$column_map  CSV with columns lab_column, standard (see config/column_map_whonet.csv)
#         data/organism_lookup.csv  reviewed overrides for names as.mo() gets wrong
# Output: cfg$processed, <cfg$reports>/cleaning_log.md, <cfg$reports>/tables/alerts.csv
#
# Antibiotic columns are recognised by name: an AMR/WHONET code (CIP), a WHONET result
# column (CIP_ND5 = CLSI disk 5 ug, CIP_EM = EUCAST MIC, CIP_NE = CLSI Etest) or the full
# drug name (Ciprofloxacin). S/I/R text is kept as the lab reported it. Raw MIC and disk values
# are interpreted with cfg$guideline, so one system decides every S/I/R. A MIC means the same
# under CLSI and EUCAST; a zone diameter does not, because the disks differ (CLSI ampicillin
# 10 ug, EUCAST 2 ug). Zones from another system's disks (CIP_ND5 under EUCAST) are therefore
# not interpreted and are counted in the log. A plain column (CIP, Ciprofloxacin) holding
# numbers stops the run: a number alone is not a MIC or a zone.

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
})

# `icr` is the D-test for inducible clindamycin resistance; a positive test makes clindamycin R.
std_cols <- c("date", "patient", "age", "gender", "ward", "specimen", "organism", "icr")
required_cols <- c("date", "patient", "organism")
# Drug code, breakpoint system (N = CLSI, E = EUCAST), method, disk content (SXT_ND1_2 = 1.25 ug).
whonet_pattern <- "^([A-Z0-9]+)_([A-Z])([DME])[0-9._]*$"

# "SENSITIVE" is a common lab synonym that as.sir() does not recognise.
clean_sir <- function(x) {
  x <- toupper(trimws(x))
  x[x == "SENSITIVE"] <- "S"
  x[grepl("[0-9]", x)] <- NA  # measurements, never S/I/R codes (as.sir() would guess)
  suppressWarnings(as.sir(x))
}

# yyyy-mm-dd or dd/mm/yyyy, whole string only. Anything else stays NA and stops the run.
iso_pattern <- "^\\d{4}-\\d{2}-\\d{2}$"
dmy_pattern <- "^\\d{2}/\\d{2}/\\d{4}$"
clean_date <- function(x) {
  x <- trimws(x)
  out <- as.Date(rep(NA_character_, length(x)))
  iso <- grepl(iso_pattern, x)
  dmy <- grepl(dmy_pattern, x)
  out[iso] <- as.Date(x[iso], format = "%Y-%m-%d")
  out[dmy] <- as.Date(x[dmy], format = "%d/%m/%Y")
  out
}

# Lookup first (case-insensitive), as.mo() for the rest. Every as.mo() match that was not
# exact is returned in `uncertain` so a person can check it and, if wrong, add it to the lookup.
clean_organism <- function(x, lookup) {
  x <- trimws(x)
  hit <- match(tolower(x), tolower(trimws(lookup$raw_name)))
  todo <- unique(x[is.na(hit) & !is.na(x)])
  mo_todo <- suppressWarnings(suppressMessages(as.mo(todo, info = FALSE)))
  unc <- suppressMessages(mo_uncertainties())
  mo <- as.character(mo_todo)[match(x, todo)]
  mo[!is.na(hit)] <- lookup$mo[hit[!is.na(hit)]]
  uncertain <- if (nrow(unc) == 0) {
    tibble(raw_name = character(), mo = character(), fullname = character())
  } else {
    tibble(raw_name = unc$original_input, mo = as.character(unc$mo), fullname = unc$fullname) |>
      filter(raw_name %in% todo)
  }
  list(mo = as.mo(mo, info = FALSE), uncertain = uncertain)
}

# Antibiotic code and test method for each column name; NA code if the name is not an
# antibiotic. Matching is exact on purpose: as.ab() guesses, and turns "Country" into
# metronidazole and "Ward" into prulifloxacin.
antibiotic_columns <- function(cols) {
  whonet <- grepl(whonet_pattern, cols)
  code <- ifelse(whonet, sub(whonet_pattern, "\\1", cols), cols)
  ab <- as.character(suppressWarnings(as.ab(code, info = FALSE)))
  ok <- !is.na(ab) & (ab == code | tolower(code) == tolower(ab_name(ab, language = NULL)))
  method <- unname(c(D = "disk", M = "mic", E = "mic")[sub(whonet_pattern, "\\3", cols)])
  system <- sub(whonet_pattern, "\\2", cols)
  system <- ifelse(system == "N", "CLSI", ifelse(system == "E", "EUCAST", system))
  tibble(column = cols, ab = ifelse(ok, ab, NA), method = ifelse(whonet, method, NA),
         system = ifelse(whonet, system, NA))
}

# S/I/R text as reported; raw MIC or disk values interpreted with `guideline`. Values with no
# breakpoint for that organism and drug stay NA.
# Any value with a digit is a measurement. as.mic()/as.disk() read signs, decimal commas and
# "=<"; combinations keep their first part (TZP "16/4" -> 16), the convention for breakpoints.
# `guideline` NA means the values cannot be read with the run's guideline; they are skipped.
interpret_column <- function(x, mo, ab, method, guideline, column = ab) {
  sir <- clean_sir(x)
  num <- !is.na(x) & grepl("[0-9]", x)
  invalid <- character()
  no_breakpoint <- 0
  if (any(num) && is.na(method)) {
    stop("Column `", column, "` holds numbers but its name does not say MIC or disk. ",
         "Rename it in the column map to WHONET style, e.g. ", ab, "_EM (MIC) or ", ab, "_ED (disk).")
  }
  if (any(num) && is.na(guideline)) {
    return(list(sir = sir, n_numeric = 0, invalid = invalid, n_no_breakpoint = 0, n_skipped = sum(num)))
  }
  if (any(num)) {
    raw_values <- trimws(sub("\\s*(mg/l|ug/ml|\u00b5g/ml|mm)$", "", x[num], ignore.case = TRUE))
    values <- suppressWarnings(if (method == "mic") as.mic(raw_values) else as.disk(raw_values))
    invalid <- unique(trimws(x[num][is.na(values)]))
    sir[num] <- suppressWarnings(suppressMessages(
      as.sir(values, mo = mo[num], ab = ab, guideline = guideline, info = FALSE)))
    no_breakpoint <- sum(!is.na(values) & is.na(sir[num]))
  }
  list(sir = sir, n_numeric = sum(num), invalid = invalid, n_no_breakpoint = no_breakpoint, n_skipped = 0)
}

# Keyed hash (HMAC-SHA256) of each patient ID. Without the secret salt the IDs cannot be
# recovered by hashing every possible hospital number; with it the lab can still find a patient.
pseudonymise <- function(id, salt) {
  u <- unique(id[!is.na(id) & id != ""])
  h <- vapply(u, function(i) substr(digest::hmac(salt, i, "sha256"), 1, 16), character(1))
  unname(h[match(id, u)])
}

read_salt <- function(path) {
  if (!file.exists(path)) {
    if (!file.exists("/dev/urandom")) {
      stop("Create ", path, " containing a long random secret (at least 32 characters).")
    }
    con <- file("/dev/urandom", "rb", raw = TRUE)
    writeLines(paste(readBin(con, "raw", 32), collapse = ""), path)
    close(con)
    message("Created ", path, ". Keep it with the data and out of git: without it, new runs ",
            "give new pseudonyms that cannot be linked to earlier ones.")
  }
  salt <- readLines(path, n = 1, warn = FALSE)
  if (nchar(salt) < 32) stop(path, " must hold at least 32 characters.")
  salt
}

clean_isolates <- function(raw, lookup, column_map = NULL, guideline = "EUCAST 2025",
                           salt = "test-salt-not-secret-0123456789ab") {
  ignored <- character()
  if (!is.null(column_map)) {
    m <- filter(column_map, lab_column %in% names(raw))
    ignored <- m$lab_column[m$standard == "ignore"]
    m <- filter(m, standard != "ignore")
    raw <- raw[setdiff(names(raw), ignored)]
    names(raw)[match(m$lab_column, names(raw))] <- m$standard
  }
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  for (col in setdiff(std_cols, names(raw))) raw[[col]] <- NA_character_

  abc <- antibiotic_columns(setdiff(names(raw), std_cols))
  if (anyNA(abc$ab)) {
    stop("Columns that are not antibiotics: ", paste(abc$column[is.na(abc$ab)], collapse = ", "),
         ". Map them in the column map, or set their standard name to `ignore`.")
  }

  organism_raw <- trimws(raw$organism)
  org <- clean_organism(organism_raw, lookup)
  mo <- org$mo
  dates <- clean_date(raw$date)
  if (anyNA(dates)) stop("Unparsed dates: ", paste(head(unique(raw$date[is.na(dates)]), 10), collapse = ", "))

  family <- sub(" .*", "", guideline)
  other_disk <- abc$method %in% "disk" & !is.na(abc$system) & abc$system != family
  abc$guideline <- ifelse(other_disk, NA, guideline)
  interpreted <- Map(function(col, ab, method, g) interpret_column(raw[[col]], mo, ab, method, g, col),
                     abc$column, abc$ab, abc$method, abc$guideline)
  sir_cols <- lapply(interpreted, `[[`, "sir")
  given <- lapply(raw[abc$column], function(r) !is.na(r) & trimws(r) != "")
  unparsed_sir <- unlist(Map(function(r, s, g) unique(trimws(r[g & is.na(s) & !grepl("[0-9]", r)])),
                             raw[abc$column], sir_cols, given))

  # Several columns can hold the same drug (CIP_ND5 and CIP_NM). The first result is kept;
  # rows where they disagree are counted.
  conflicts <- 0
  merged <- list()
  for (i in seq_along(sir_cols)) {
    ab <- abc$ab[i]
    s <- sir_cols[[i]]
    if (is.null(merged[[ab]])) {
      merged[[ab]] <- s
    } else {
      conflicts <- conflicts + sum(!is.na(merged[[ab]]) & !is.na(s) & merged[[ab]] != s)
      merged[[ab]][is.na(merged[[ab]])] <- s[is.na(merged[[ab]])]
    }
  }

  # A low clindamycin MIC misses inducible resistance in staphylococci; labs report R after a
  # positive D-test, and so must the cleaned data.
  icr_positive <- toupper(trimws(raw$icr)) %in% c("+", "POS", "POSITIVE", "YES", "Y", "TRUE", "1")
  icr_changed <- 0
  if (!is.null(merged$CLI)) {
    icr_changed <- sum(icr_positive & (is.na(merged$CLI) | merged$CLI != "R"))
    merged$CLI[icr_positive] <- "R"
  }

  clean <- tibble(
    date = dates,
    patient = pseudonymise(trimws(raw$patient), salt),
    age = suppressWarnings(as.integer(raw$age)),
    gender = toupper(trimws(raw$gender)),
    ward = trimws(raw$ward),
    specimen = trimws(raw$specimen),
    mo = mo,
    !!!merged
  )
  n_before <- nrow(clean)
  clean <- distinct(clean)

  log <- list(
    rows_in = nrow(raw),
    rows_out = nrow(clean),
    exact_duplicates_removed = n_before - nrow(clean),
    ignored_columns = ignored,
    antibiotic_columns = abc,
    guideline = guideline,
    numeric_values = sum(vapply(interpreted, `[[`, numeric(1), "n_numeric")),
    skipped = tibble(column = abc$column, n = vapply(interpreted, `[[`, numeric(1), "n_skipped")) |>
      filter(n > 0),
    no_breakpoint = tibble(column = abc$column,
                           n = vapply(interpreted, `[[`, numeric(1), "n_no_breakpoint")) |>
      filter(n > 0) |> arrange(desc(n)),
    invalid_values = unique(unlist(lapply(interpreted, `[[`, "invalid"))),
    same_drug_conflicts = conflicts,
    icr_changed = icr_changed,
    organism_spellings = length(unique(organism_raw)),
    organisms_matched = length(unique(mo[!is.na(mo) & mo != "UNKNOWN"])),
    overridden_by_lookup = sum(tolower(organism_raw) %in% tolower(trimws(lookup$raw_name))),
    unmatched_organisms = unique(organism_raw[is.na(mo) | mo == "UNKNOWN"]),
    uncertain_organisms = org$uncertain,
    sir_values_changed = sum(unlist(Map(function(r, s, g) sum(g & !grepl("[0-9]", r) &
                                                              r != as.character(s), na.rm = TRUE),
                                        raw[abc$column], sir_cols, given))),
    unparsed_sir = unparsed_sir,
    dmy_dates = sum(grepl(dmy_pattern, trimws(raw$date)))
  )
  list(data = clean, log = log)
}

# Isolates infection control needs to hear about. One row per isolate and alert.
# Intrinsic resistance never raises an alert: E. gallinarum and E. casseliflavus carry vanC and
# are vancomycin resistant by nature, which is not VRE. CRE follows the CDC definition: resistant to any carbapenem, except imipenem alone for
# Proteus, Providencia and Morganella, which have naturally raised imipenem MICs.
flag_alerts <- function(data) {
  any_r <- function(codes, rows = TRUE) {
    codes <- intersect(codes, names(data))
    hit <- Reduce(`|`, lapply(codes, function(ab) {
      !is.na(data[[ab]]) & data[[ab]] == "R" & !mo_is_intrinsic_resistant(data$mo, ab)
    }), FALSE)
    hit & rows
  }
  drugs_r <- function(codes) {
    codes <- intersect(codes, names(data))
    apply(as.data.frame(lapply(data[codes], function(x) !is.na(x) & x == "R")), 1,
          function(r) paste(codes[r], collapse = ", "))
  }
  genus <- mo_genus(data$mo, language = NULL)
  species <- mo_species(data$mo, language = NULL)
  carbapenems <- c("IPM", "MEM", "ETP", "DOR")
  low_ipm <- genus %in% c("Proteus", "Providencia", "Morganella")
  rules <- list(
    MRSA = list(genus %in% "Staphylococcus" & species %in% "aureus", c("OXA", "FOX", "FLC", "MET")),
    VRE = list(genus %in% "Enterococcus", "VAN"),
    CRE = list(mo_order(data$mo, language = NULL) %in% "Enterobacterales", carbapenems),
    `CR Acinetobacter baumannii` = list(genus %in% "Acinetobacter" & grepl("^baumannii", species),
                                        carbapenems),
    `CR Pseudomonas aeruginosa` = list(genus %in% "Pseudomonas" & species %in% "aeruginosa", carbapenems)
  )
  empty <- tibble(alert = character(), date = as.Date(character()), patient = character(),
                  ward = character(), specimen = character(), organism = character(),
                  resistant_to = character())
  bind_rows(empty, lapply(names(rules), function(alert) {
    rows <- rules[[alert]][[1]]
    codes <- rules[[alert]][[2]]
    hit <- if (alert == "CRE") any_r(setdiff(codes, "IPM"), rows) | any_r("IPM", rows & !low_ipm) else any_r(codes, rows)
    if (!any(hit)) return(NULL)
    tibble(alert = alert, date = data$date[hit], patient = data$patient[hit], ward = data$ward[hit],
           specimen = data$specimen[hit], organism = mo_name(data$mo[hit], language = NULL),
           resistant_to = drugs_r(codes)[hit])
  })) |> arrange(alert, date)
}

write_cleaning_log <- function(log, path) {
  fmt <- function(x) if (length(x) == 0) "none" else paste0("`", x, "`", collapse = ", ")
  abc <- log$antibiotic_columns
  lines <- c(
    "# Cleaning log",
    "",
    "Generated by `R/01_clean.R`. Do not edit by hand.",
    "",
    sprintf("- Rows in: %d", log$rows_in),
    sprintf("- Exact duplicate rows removed: %d", log$exact_duplicates_removed),
    sprintf("- Rows out: %d", log$rows_out),
    "- Patient IDs replaced by keyed hashes (HMAC-SHA256 with the secret salt file)",
    sprintf("- Columns dropped (set to `ignore` in the column map): %s", fmt(log$ignored_columns)),
    sprintf("- Dates converted from dd/mm/yyyy: %d", log$dmy_dates),
    sprintf("- Distinct organism spellings: %d, standardised with `as.mo()` to %d organisms",
            log$organism_spellings, log$organisms_matched),
    sprintf("- Names set from `data/organism_lookup.csv` (reviewed overrides): %d rows",
            log$overridden_by_lookup),
    sprintf("- Names that became `UNKNOWN`: %s", fmt(log$unmatched_organisms)),
    sprintf("- S/I/R values rewritten (case, spaces, words to codes): %d", log$sir_values_changed),
    sprintf("- Unparsed S/I/R values (set to NA): %s", fmt(log$unparsed_sir)),
    sprintf("- Raw MIC and disk values interpreted with %s breakpoints: %d", log$guideline, log$numeric_values),
    sprintf("- Disk zones from another system's disks, not interpreted (set to NA): %d in %s",
            sum(log$skipped$n), fmt(log$skipped$column)),
    sprintf("- Of those, values with no breakpoint for that organism and drug (set to NA): %d",
            sum(log$no_breakpoint$n)),
    sprintf("- Numbers that are not valid MIC or disk values (set to NA): %s", fmt(log$invalid_values)),
    sprintf("- Clindamycin set to R after a positive D-test (inducible resistance): %d", log$icr_changed),
    sprintf("- Isolates where two columns for the same drug disagree (first column kept): %d",
            log$same_drug_conflicts),
    "",
    "Exact duplicates are rows identical in every column. Without an isolate ID they cannot be",
    "told apart from repeat tests, and `first_isolate()` would drop them later anyway.",
    "",
    "## Antibiotic columns",
    "",
    "| Column | Drug | Numbers read as |",
    "|---|---|---|",
    sprintf("| `%s` | %s (`%s`) | %s |", abc$column, ab_name(abc$ab, language = NULL), abc$ab,
            ifelse(is.na(abc$method), "S/I/R as reported",
                   ifelse(is.na(abc$guideline), paste(abc$method, "from another system, not read"),
                          paste(abc$method, "with", abc$guideline)))),
    if (nrow(log$no_breakpoint) > 0) c(
      "",
      "## Values without a breakpoint",
      "",
      "The guideline has no breakpoint for some organism and drug combinations in these columns.",
      "Those values are NA in the cleaned data, not susceptible.",
      "",
      "| Column | Values |",
      "|---|---|",
      sprintf("| `%s` | %d |", log$no_breakpoint$column, log$no_breakpoint$n)
    ),
    "",
    "## Organism names matched by approximation",
    "",
    "`as.mo()` matched these names without an exact hit. Check each one. If a match is wrong,",
    "add the raw name and the correct code to `data/organism_lookup.csv` and run again.",
    "Abbreviations such as `S. hominis` are ambiguous (Staphylococcus or Streptococcus?).",
    "",
    "| Raw name | Matched code | Matched name |",
    "|---|---|---|",
    sprintf("| %s | `%s` | %s |", log$uncertain_organisms$raw_name,
            log$uncertain_organisms$mo, log$uncertain_organisms$fullname)
  )
  writeLines(lines, path)
}

if (!isTRUE(getOption("amrwf.skip_main"))) {
  if (!exists("cfg")) source("config/demo.R")
  if (isTRUE(cfg$demo) && cfg$input != "data/raw/isolates_messy.csv") {
    stop("cfg$demo is TRUE but the input is not the demo file. Set demo = FALSE and give the ",
         "export its own `processed` and `reports` paths (see config/demo.R).")
  }
  raw <- read_csv(cfg$input, col_types = cols(.default = col_character()), na = "")
  lookup <- read_csv("data/organism_lookup.csv", col_types = cols(.default = col_character()))
  column_map <- if (is.null(cfg$column_map)) NULL else
    read_csv(cfg$column_map, col_types = cols(.default = col_character()))
  result <- clean_isolates(raw, lookup, column_map, cfg$guideline, read_salt(cfg$salt_file))
  dir.create(file.path(cfg$reports, "tables"), showWarnings = FALSE, recursive = TRUE)
  write_csv(flag_alerts(result$data), file.path(cfg$reports, "tables/alerts.csv"), na = "")
  out <- result$data |> mutate(mo = as.character(mo), across(where(is.sir), as.character))
  write_csv(out, cfg$processed, na = "")
  write_cleaning_log(result$log, file.path(cfg$reports, "cleaning_log.md"))
  message("Wrote ", cfg$processed, " (", nrow(out), " rows), cleaning log and alerts")
}
