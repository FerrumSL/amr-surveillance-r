# Workflow guide

A short guide to the scripts in `R/`. Run everything with `Rscript R/run_all.R [config file]`
from the project root after `renv::restore()`. All settings live in the config file
(`config/demo.R` by default); the scripts should not need editing.

## Purpose

Turn one laboratory export of isolate test results into standard surveillance outputs:
a cleaned dataset with a log, antibiograms, resistance by period with confidence intervals,
and a short written summary. The scripts use the
[AMR R package](https://amr-for-r.org) (Berends et al.) for organism names, S/I/R
interpretation, first-isolate selection, antibiograms and intrinsic resistance.

## Data assumptions

- One row per isolate. Required columns `date`, `patient`, `organism`; optional `age`,
  `gender`, `ward`, `specimen` and `icr` (D-test for inducible clindamycin resistance). Other
  names are translated by the column map (`lab_column,standard`; `ignore` drops a column).
- One column per antibiotic, named by its code (`CIP`), its WHONET result name (`CIP_ND5`,
  `CIP_EM`, `CIP_NE`: N/E = CLSI/EUCAST, D = disk, M/E = MIC/Etest) or its full name
  (`Ciprofloxacin`). Matching is exact; anything else stops the run.
- Results are S/I/R text or raw values. Text accepts case and spacing differences and the
  words susceptible, sensitive, intermediate and resistant, and is kept as the lab reported
  it. MIC and disk values are interpreted with `cfg$guideline` (e.g. `EUCAST 2025`). Values
  with no breakpoint stay NA and are counted in the cleaning log. When two columns hold the
  same drug (`CIP_ND5` and `CIP_EE`), the first result is kept and disagreements are counted.
- Dates are `yyyy-mm-dd` or `dd/mm/yyyy`. Any other format stops the run, because a guess
  between day-first and month-first could silently swap months.
- The demo file `data/raw/isolates_messy.csv` is the package's `example_isolates` data
  (fictitious) with mess added on purpose by `R/00_fetch.R`.

## Steps and why

1. **`R/01_clean.R`: clean, pseudonymise, flag.**
   - Patient IDs become HMAC-SHA256 hashes keyed with `salt.txt` (kept out of git). Only
     standard and antibiotic columns are written, so names and comments never reach
     `data/processed/`.
   - Infection-control alerts (`tables/alerts.csv`): MRSA (oxacillin, cefoxitin or
     flucloxacillin R), VRE, CRE (CDC definition), carbapenem-resistant *A. baumannii* and
     *P. aeruginosa*. All isolates, not first isolates; intrinsic resistance never alerts.
   - Organism names go through `data/organism_lookup.csv` first, then `as.mo()`.
     Abbreviations are ambiguous: `as.mo("S. hominis")` returns *Streptococcus* hominis, but
     the lab meant *Staphylococcus* hominis. Every approximate match is listed in
     the cleaning log (`reports/<run>/cleaning_log.md`) for a person to check; wrong ones go into the lookup.
   - Results go through `as.sir()`. Values it cannot read become NA and are listed in the log.
   - Exact duplicate rows are dropped and counted.
   - The test `tests/testthat/test_clean.R` checks that cleaning the messy file gives back the
     package's original data value for value.
2. **`R/02_analyse.R`: analyse.** Parameters are at the top of the file.
   - **First isolates.** Only the first isolate of each species per patient per 365 days is
     kept (`method = "episode-based"`), following CLSI M39. Without this, a patient sampled
     ten times during one infection would count ten times, and resistant infections, which
     last longer and are sampled more often, would inflate resistance.
   - **Guideline.** EUCAST counts I ("susceptible, increased exposure") as susceptible; CLSI
     counts it as resistant. The choice is set once (`guideline`) and used everywhere. If the
     export has raw MIC or zone diameters, interpret them with `as.sir()` and a stated
     guideline version, and never mix versions across years without saying so.
   - **Antibiograms** (`antibiogram()`) for `antibiogram_years`, split into Gram-negative and
     Gram-positive panels, also printed to `antibiogram.pdf` (A4, CLSI M39 layout: organism
     with number of first isolates, % susceptible per drug).
     Intrinsic resistance (EUCAST Expected Resistant Phenotypes) is blanked, because it is
     biology, not a surveillance signal.
   - **Resistance by period** with exact (Clopper-Pearson) 95% confidence intervals. No
     percentage is shown for fewer than 30 tested isolates.
   - **Trend test**: logistic regression of resistant (yes/no) on calendar year, one row per
     first isolate, giving an odds ratio per year with a profile-likelihood CI and a
     likelihood-ratio p-value.
3. **`R/03_report.R`: summary.** Writes `summary.md`. Numbers are inserted by the script.
   For the demo, the interpretation sentences are hand-written and each claim they make is
   checked with `stopifnot()`. For other data the summary holds the numbers and alert counts
   and leaves the interpretation to a person.
4. **`R/04_validate.R`: interpretation check** (NCBI data only). Re-interprets real MIC and
   disk values and compares them with the S/I/R the submitting labs reported, using the CLSI
   M52 error classes.

## Running it on a new export

1. Copy `config/demo.R` to `config/my_lab.R`. Set `demo = FALSE`, `input`, `guideline`,
   `processed`, `reports`, the years and the antibiotic panels.
2. If the export's column names differ, write a column map (start from
   `config/column_map_whonet.csv` for WHONET) and set `column_map`.
3. Run `Rscript R/01_clean.R` with `cfg` loaded, or the whole run, and read the cleaning log.
   Fix wrong organism matches in `data/organism_lookup.csv` and run again until it is clean.
4. Run `Rscript R/run_all.R config/my_lab.R`, then write the interpretation by hand.

## Known pitfalls

- **Small numbers.** A 100% resistance figure from 4 isolates says almost nothing. Always show
  n and the confidence interval; the minimum of 30 is a floor, not a guarantee.
- **Intrinsic resistance** must not be reported as acquired resistance.
- **Changing guidelines.** EUCAST and CLSI breakpoints change most years. A trend across a
  breakpoint change can be an artefact. Record the version used for each year.
- **Changing lab practice.** In the demo data, coagulase-negative staphylococci are reported
  as "CoNS" until 2012 and as species from 2013. Trends by species across that change are
  not real.
- **Multiple testing.** Testing many pairs will produce some p < 0.05 by chance.
- **Testing bias.** Second-line drugs are often tested only on isolates already resistant to
  first-line drugs, which inflates their resistance. Compare `n_tested` with the number of
  isolates.
- **Combination MICs** (`16/4`) are read by their first part, the convention breakpoints use.
- **Drugs without an S category** (colistin under CLSI) print as n/a in the antibiogram, not 0%.
- **Inducible clindamycin resistance** is invisible in a MIC. Map the D-test column to `icr`.
- **Selected isolates.** Collections like NCBI BioSample are chosen by their submitters and
  over-represent resistance. Use them to test the pipeline, not to estimate rates.
