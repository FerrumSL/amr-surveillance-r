# AMR surveillance workflow in R

[![run-all](https://github.com/FerrumSL/amr-surveillance-r/actions/workflows/run-all.yml/badge.svg)](https://github.com/FerrumSL/amr-surveillance-r/actions/workflows/run-all.yml)

R scripts that take a lab export of antibiotic test results and turn it into the usual
surveillance outputs: a cleaned, pseudonymised dataset, a one-page printable antibiogram,
resistance over time with confidence intervals, an infection-control alert list (MRSA, VRE,
CRE, carbapenem-resistant *A. baumannii* and *P. aeruginosa*) and a short written summary.
The heavy lifting is done by the [AMR package](https://amr-for-r.org).

It reads exports with any column names (through a column map), WHONET exports (a ready
map is included), and S/I/R results or raw MIC and disk-zone values, which it interprets with
one stated EUCAST or CLSI version.

Reading every year with that one version keeps its trends free of breakpoint revisions, but it
cannot say how much the antibiograms a lab published at the time were moved by them.
[breakpoint-drift](https://github.com/FerrumSL/breakpoint-drift) does that split.

## Example: public NCBI isolates

![Cumulative antibiogram, 1998-2025, of 2,004 human isolates from NCBI BioSample interpreted with CLSI 2025: percentage susceptible for A. baumannii, E. coli, K. pneumoniae, P. aeruginosa, E. faecium and S. aureus by drug grouped by class, shaded in six blue bands from 90% or more (lightest) to under 20% (darkest); IR marks intrinsic resistance](reports/ncbi/figures/antibiogram.png)

This page comes from 2,138 real human isolates in
[NCBI BioSample](https://www.ncbi.nlm.nih.gov/biosample/) antibiograms (*E. coli*,
*K. pneumoniae*, *S. aureus*, *P. aeruginosa*, *A. baumannii*, *E. faecium*), with raw MIC and disk
values and the S/I/R each submitting lab reported. The pipeline reads the raw values with CLSI 2025
and checks its reading against the labs' own calls (see below). These isolates were picked by
their submitters, often because they were resistant, so their resistance rates say nothing about
any population and must not guide therapy.

Outputs in [`reports/ncbi/`](reports/ncbi/):
[`antibiogram.pdf`](reports/ncbi/antibiogram.pdf) (one A4 page),
[`summary.md`](reports/ncbi/summary.md), [`cleaning_log.md`](reports/ncbi/cleaning_log.md),
[`validation.md`](reports/ncbi/validation.md) (the check against the labs' calls), and
`tables/` with the antibiograms, resistance by period, trend tests and `alerts.csv` as CSV.

## Synthetic test data

The default run (`config/demo.R`) uses `example_isolates`, synthetic data that ships with the
AMR package: 2,000 made-up isolates from 2002 to 2017. It covers what the NCBI data cannot: real
patients with repeat isolates (so first-isolate selection matters), wards for the alert list,
and a messy export. `R/00_fetch.R` messes it up on purpose (organism names spelled several ways,
results written as "s", "Resistant" or "SENSITIVE", some dates in dd/mm/yyyy), and a test checks
that cleaning recovers the original data exactly. It is also the only run with a hand-written
interpretation, which `R/03_report.R` checks claim by claim against the numbers. None of its
numbers describe a real hospital or country.

Outputs in [`reports/demo/`](reports/demo/): [`antibiogram.pdf`](reports/demo/antibiogram.pdf),
[`summary.md`](reports/demo/summary.md), [`cleaning_log.md`](reports/demo/cleaning_log.md) and
`tables/`. Its trend chart is the clearer one, because the NCBI isolates are too few per period
for a trend:

![Line charts of the percentage of E. coli first isolates resistant to amoxicillin/clavulanic acid, ciprofloxacin and trimethoprim/sulfamethoxazole in four 4-year periods from 2002 to 2017, with 95% confidence intervals and the number tested per period (synthetic test data)](reports/demo/figures/resistance_trend.png)

## How the antibiogram is drawn

The antibiogram uses one blue scale in six bands whose lightness differs by at least 12 L\*,
so it reads the same printed in black and white. Text in each cell is black or white,
whichever reaches at least 4.5:1 contrast. `R/plot_style.R` holds the palette and the checks,
and the run stops if either fails. Grayscale proofs are written to `reports/_grayscale_check/`
(not committed). Cells marked IR are intrinsic resistance, a dash means fewer than 30
isolates tested, and a blank means not tested.

## How to run it

You need R (tested on 4.6.1). From the project folder:

```sh
Rscript -e 'renv::restore()'
Rscript R/run_all.R
```

This cleans the data, runs the analysis, writes the reports and then runs the tests.
`Rscript R/run_all.R config/ncbi.R` does the same for the real NCBI isolates.

### On your own lab's export

1. Copy `config/demo.R` to `config/my_lab.R`. Set `demo = FALSE`, `name`, `input`,
   `guideline`, the years, and new `processed` and `reports` paths so the demo outputs are not
   overwritten. The run stops if `demo = TRUE` is left on with another input.
2. If the column names differ from `date`, `patient`, `organism` (and optionally `age`,
   `gender`, `ward`, `specimen`), point `column_map` at a CSV with columns `lab_column` and
   `standard`. Set `standard` to `ignore` for columns to drop (names, comments). For WHONET
   exports use `config/column_map_whonet.csv`.
3. Antibiotic columns are found by name: a code (`CIP`), a WHONET column (`CIP_ND5` = CLSI
   disk, `CIP_EM` = EUCAST MIC) or the drug name (`Ciprofloxacin`). Every value is read with
   `guideline`. Disk zones from the other system's disks (`_ND` under EUCAST, `_ED` under
   CLSI) are not read, because the disk contents differ; the cleaning log counts them. A
   column of numbers without MIC or disk in its name stops the run; rename it in the column
   map (`CIP_EM`).
4. Run `Rscript R/run_all.R config/my_lab.R` and read the `cleaning_log.md` in your `reports`
   folder first.

Patient IDs are replaced with keyed hashes (HMAC-SHA256) before anything is saved. The key is
in `salt.txt`, created on the first run and kept out of git. Keep it: the same key gives the
same pseudonym on every run, which is what lets first isolates be counted across years, and
it lets the lab look a patient up again. Without it the IDs cannot be recovered by hashing
every possible hospital number. Columns not mapped to a standard name are never written out.

[`docs/WORKFLOW.md`](docs/WORKFLOW.md) walks through each step, why it's done that way, and
how to point the scripts at a different export.

## What's in here

```
config/demo.R          settings for the demo run (copy it for a new export)
config/ncbi.R          settings for the real NCBI isolates
config/column_map_*.csv  export column names -> standard names (WHONET, NCBI)
R/00_fetch.R           exports the example data and adds the deliberate mess
R/00_fetch_ncbi.R      downloads the NCBI antibiograms and builds a WHONET-style export
R/01_clean.R           column map, organism names, S/I/R and MIC/disk interpretation,
                       dates, duplicates, pseudonymisation, infection-control alerts
R/02_analyse.R         first isolates, antibiograms (CSV, PNG, PDF), resistance by period, trend test
R/03_report.R          writes summary.md; for the demo, checks its claims against the numbers
R/plot_style.R         palette, theme, and the L* and contrast checks for print
R/04_validate.R        compares the pipeline's MIC/disk interpretation with the labs' calls
R/run_all.R            runs all of the above, then the tests
data/raw/              the inputs (demo, NCBI snapshot)
data/organism_lookup.csv  manual fixes for names the AMR package matches wrongly
data/processed/        the cleaned, pseudonymised isolates
tests/testthat/        tests for the cleaning, interpretation, WHONET import and alerts
```

## A few things worth knowing

- Abbreviated names are a trap. `as.mo("S. hominis")` returns *Streptococcus* hominis, but a
  lab writing that almost always means *Staphylococcus* hominis. I found four of these and
  fixed them in `data/organism_lookup.csv`. The cleaning log lists every name that was matched
  approximately, so they can be checked by hand.
- The test suite cleans the messy file and checks that the result matches the package's
  original data exactly, value for value.
- Organisms that are naturally resistant to a drug (for example *K. pneumoniae* and
  amoxicillin) are left blank in the antibiograms, since that isn't acquired resistance.
- I used `example_isolates` rather than the package's `WHONET` example for the demo, because
  `WHONET` only has 500 isolates spread over two countries and 16 years, which is too few for
  trends. The tests use `WHONET` to check the WHONET import.

## Checks built in

Each of these is a common way an AMR pipeline goes wrong on real exports. The real NCBI data
were used to test them, and most have a unit test.

- Antibiotic columns are matched exactly. `as.ab()` alone would read a `Country` column as
  metronidazole and `Ward` as prulifloxacin.
- Disk zones are read only with breakpoints for the same disks: CLSI and EUCAST disks differ in
  content, so zones from the other system are counted and left out (595 in the NCBI data).
- Combination MICs (`16/4`, `0.5/9.5`) are read by their first part; values that are not valid
  MICs or zones are listed in the cleaning log, never dropped quietly.
- Values with no breakpoint for that organism and drug stay blank and are counted per column;
  they are never read as susceptible (about 3,100 NCBI values, mostly veterinary drugs).
- Not-interpretable (NI) results are kept out of denominators.
- Intrinsic resistance never raises an alert (no VRE alert for vanC *E. casseliflavus*), and CRE
  ignores imipenem alone for *Proteus*, *Providencia* and *Morganella*.
- A positive D-test sets clindamycin to R. When there is no D-test result, the antibiogram marks
  staphylococcal clindamycin with ‡ and says how many erythromycin-R, clindamycin-S isolates
  could hide inducible resistance (120 in the NCBI data).
- Drugs without an S category in the guideline (colistin under CLSI) print n/a, not 0%.

[`reports/ncbi/validation.md`](reports/ncbi/validation.md) compares the pipeline's CLSI 2025
reading of 15,328 real MIC and disk values with the labs' own calls: 94.6% agree. With the CLSI
version of each isolate's collection year (possible for 10,603 tests), major errors (pipeline R,
lab S) fall from 52 to 18, so many of the disagreements come from breakpoint changes. The largest
gap is *S. aureus* and clindamycin (127 isolates the pipeline calls S were reported R), the
inducible-resistance pattern above.

## Known limitations

- **Selected isolates.** NCBI submitters chose what to sequence, often resistant strains, so
  the NCBI rates describe no population. Use them to test the pipeline, not to guide therapy.
- **Completed dates.** 867 of the 2,004 NCBI isolates have only a month (209) or a year (658);
  their date is set to the 15th or 1 July, recorded in `Date precision`. This is safe here because
  each isolate is its own patient and every analysis works by year. Lab exports are never
  completed: an inexact date stops the run.
- **No D-test in the NCBI data.** The *S. aureus* clindamycin figure there may be too high
  (marked ‡ on the page).
- **Intrinsic resistance follows EUCAST** (Expected Resistant Phenotypes v1.2) even on the CLSI
  page, because the AMR package has no copy of CLSI M100 Appendix B. The two lists mostly agree.
- **One page per run.** An antibiogram with more than about 12 organism rows needs a second page,
  which is not built yet.

## Licence and credits

The code is released under the GNU General Public License v2 (see [`LICENSE`](LICENSE)).

The demo data (`data/raw/isolates_messy.csv`, `data/processed/isolates_clean.csv`) is derived
from the `example_isolates` dataset in the AMR package, which is also GPL-2. Intrinsic resistance rules come from EUCAST Expected Resistant Phenotypes
v1.2 (2023), as bundled in the AMR package.

The NCBI isolates in `data/raw/ncbi_*.csv` come from NCBI BioSample antibiograms,
downloaded on 2026-10-21 by `R/00_fetch_ncbi.R`. NCBI places no restrictions on their use,
although individual submitters may; each isolate's accession is kept in the raw files.

If you use the AMR package, please cite it:

> Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C (2022). AMR: An R Package
> for Working with Antimicrobial Resistance Data. *Journal of Statistical Software*, 104(3),
> 1-31. [doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)
