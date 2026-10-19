# Workflow guide

A short guide to the four scripts in `R/`. Run everything with `Rscript R/run_all.R` from the
project root after `renv::restore()`.

## Purpose

Turn one laboratory export of isolate test results into standard surveillance outputs:
a cleaned dataset with a log, antibiograms, resistance by period with confidence intervals,
and a short written summary. The scripts use the
[AMR R package](https://amr-for-r.org) (Berends et al.) for organism names, S/I/R
interpretation, first-isolate selection, antibiograms and intrinsic resistance.

## Data assumptions

- One row per isolate. Columns `date`, `patient`, `age`, `gender`, `ward`, `organism`, then
  one column per antibiotic named by its code (`CIP`, `AMC`, ...; `AMR::as.ab()` must
  recognise each one).
- Results are S/I/R text. The cleaner accepts case and spacing differences and the words
  susceptible, sensitive, intermediate and resistant.
- Dates are `yyyy-mm-dd` or `dd/mm/yyyy`. Any other format stops the run, because a guess
  between day-first and month-first could silently swap months.
- The demo file `data/raw/isolates_messy.csv` is the package's `example_isolates` data
  (fictitious) with mess added on purpose by `R/00_fetch.R`.

## Steps and why

1. **`R/01_clean.R`: clean and structure.**
   - Organism names go through `data/organism_lookup.csv` first, then `as.mo()`.
     Abbreviations are ambiguous: `as.mo("S. hominis")` returns *Streptococcus* hominis, but
     the lab meant *Staphylococcus* hominis. Every approximate match is listed in
     `reports/cleaning_log.md` for a person to check; wrong ones go into the lookup.
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
   - **Antibiograms** (`antibiogram()`), split into Gram-negative and Gram-positive panels.
     Intrinsic resistance (EUCAST Expected Resistant Phenotypes) is blanked, because it is
     biology, not a surveillance signal.
   - **Resistance by period** with exact (Clopper-Pearson) 95% confidence intervals. No
     percentage is shown for fewer than 30 tested isolates.
   - **Trend test**: logistic regression of resistant (yes/no) on calendar year, one row per
     first isolate, giving an odds ratio per year with a profile-likelihood CI and a
     likelihood-ratio p-value.
3. **`R/03_report.R`: summary.** Writes `reports/summary.md`. Numbers are inserted by the
   script. The interpretation sentences are hand-written, and each claim they make is
   checked with `stopifnot()`, so new data that no longer supports a sentence stops the run.

## Running it on a new country's export

1. Map the export's columns to the names above (WHONET exports: `AMR::as.ab()` reads codes
   such as `CIP_ND5`; rename them to the plain code).
2. Save it as `data/raw/isolates_messy.csv` (or change the path in `R/01_clean.R`).
3. Run `Rscript R/01_clean.R` and read `reports/cleaning_log.md`. Fix wrong organism matches
   in `data/organism_lookup.csv` and run again until the log is clean.
4. Edit the parameters at the top of `R/02_analyse.R`: years, period length, organism and
   antibiotic pairs, guideline.
5. Run `Rscript R/run_all.R`. If `R/03_report.R` stops on a claim, rewrite the interpretation.

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
