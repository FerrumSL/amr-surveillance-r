# Download real antimicrobial susceptibility test (AST) results from NCBI BioSample.
#
# Isolates from human hosts that carry an NCBI antibiogram table: the MIC or disk-zone
# measurement, the submitting lab's own S/I/R call and the testing standard it used.
# Writes two files:
#   data/raw/ncbi_ast_long.csv.gz  one row per isolate and test, as published (the snapshot)
#   data/raw/ncbi_export.csv    the same tests as a lab-style wide export, one row per isolate,
#                               antibiotic columns in WHONET style (CIP_NM = CLSI MIC)
# The pipeline reads the wide file; R/04_validate.R compares its interpretation with the
# lab calls in the long file.
#
# These are not surveillance data: submitters chose which isolates to sequence (outbreaks,
# resistant strains, reference collections), so resistance rates from them are biased upwards.
# NCBI places no restrictions on reuse; individual submitters may.

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
  library(xml2)
})

species <- c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus",
             "Pseudomonas aeruginosa", "Acinetobacter baumannii", "Enterococcus faecium")
per_species <- 400
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

search_ids <- function(organism) {
  term <- sprintf('antibiogram[filter] AND "%s"[orgn] AND "homo sapiens"[host]', organism)
  url <- paste0(eutils, "esearch.fcgi?db=biosample&retmax=", per_species,
                "&term=", URLencode(term, reserved = TRUE))
  Sys.sleep(0.4)  # NCBI allows 3 requests per second without an API key
  xml_text(xml_find_all(read_xml(url), "//IdList/Id"))
}

fetch_samples <- function(ids) {
  chunks <- split(ids, ceiling(seq_along(ids) / 200))
  unlist(lapply(chunks, function(chunk) {
    Sys.sleep(0.4)
    url <- paste0(eutils, "efetch.fcgi?db=biosample&id=", paste(chunk, collapse = ","))
    xml_find_all(read_xml(url), "//BioSample")
  }), recursive = FALSE)
}

attr_value <- function(sample, name) {
  x <- xml_find_first(sample, sprintf(".//Attribute[@harmonized_name='%s']", name))
  if (inherits(x, "xml_missing")) NA_character_ else xml_text(x)
}

parse_sample <- function(sample) {
  table <- xml_find_first(sample, ".//Table[starts-with(@class, 'Antibiogram')]")
  header <- xml_text(xml_find_all(table, "./Header/Cell"))
  rows <- xml_find_all(table, "./Body/Row")
  if (length(rows) == 0) return(NULL)
  cells <- do.call(rbind, lapply(rows, function(r) xml_text(xml_find_all(r, "./Cell"))))
  # Submitters use different versions of the table, so columns are looked up by name.
  col <- function(name) {
    i <- match(name, header)
    if (is.na(i)) rep(NA_character_, nrow(cells)) else cells[, i]
  }
  tibble(
    accession = xml_attr(sample, "accession"),
    organism = xml_text(xml_find_first(sample, ".//Organism/OrganismName")),
    collection_date = attr_value(sample, "collection_date"),
    isolation_source = attr_value(sample, "isolation_source"),
    antibiotic = col("Antibiotic"),
    phenotype = col("Resistance phenotype"),
    sign = col("Measurement sign"),
    measurement = col("Measurement"),
    units = col("Measurement units"),
    method = col("Laboratory typing method"),
    standard = col("Testing standard")
  )
}

if (!file.exists("data/raw/ncbi_ast_long.csv.gz")) {
  samples <- unlist(lapply(species, function(s) fetch_samples(search_ids(s))), recursive = FALSE)
  long <- bind_rows(lapply(samples, parse_sample))
  write_csv(long, "data/raw/ncbi_ast_long.csv.gz", na = "")
  message("Wrote data/raw/ncbi_ast_long.csv.gz (", nrow(long), " tests, ",
          n_distinct(long$accession), " isolates)")
}

# ---- Wide, WHONET-style export -----------------------------------------------------------------
long <- read_csv("data/raw/ncbi_ast_long.csv.gz", col_types = cols(.default = col_character()), na = "")

tests <- long |>
  mutate(
    ab = as.character(suppressWarnings(as.ab(antibiotic, info = FALSE))),
    # MIC in mg/L (= ug/mL) or zone diameter in mm; anything else is left out.
    kind = case_when(units == "mm" ~ "D", tolower(units) %in% c("mg/l", "ug/ml") ~ "M"),
    guideline = case_when(grepl("^CLSI", standard) ~ "N", grepl("^EUCAST", standard) ~ "E"),
    value = paste0(ifelse(sign %in% c("==", "="), "", sign), measurement),
    column = paste0(ab, "_", guideline, kind)
  ) |>
  filter(!is.na(ab), !is.na(kind), !is.na(guideline), !is.na(measurement))

isolates <- distinct(long, accession, organism, collection_date, isolation_source) |>
  filter(grepl("^\\d{4}-\\d{2}-\\d{2}$", collection_date))

wide <- tests |>
  filter(accession %in% isolates$accession) |>
  distinct(accession, column, .keep_all = TRUE) |>  # a few isolates repeat the same test
  select(accession, column, value) |>
  as.data.frame() |>
  reshape(idvar = "accession", timevar = "column", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
wide <- isolates |>
  inner_join(wide, by = "accession") |>
  rename(`Isolate ID` = accession, Organism = organism, `Collection date` = collection_date,
         `Isolation source` = isolation_source)
wide <- wide[c(names(wide)[1:4], sort(names(wide)[-(1:4)]))]

write_csv(wide, "data/raw/ncbi_export.csv", na = "")
message("Wrote data/raw/ncbi_export.csv (", nrow(wide), " isolates with a full collection date; ",
        n_distinct(long$accession) - nrow(wide), " left out without one)")
