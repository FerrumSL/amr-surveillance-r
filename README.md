# AMR surveillance workflow in R

[![run-all](https://github.com/FerrumSL/amr-surveillance-r/actions/workflows/run-all.yml/badge.svg)](https://github.com/FerrumSL/amr-surveillance-r/actions/workflows/run-all.yml)

R scripts that take a lab export of antibiotic test results and turn it into the usual
surveillance outputs: a cleaned dataset, antibiograms, resistance over time with confidence
intervals, and a short written summary. The heavy lifting is done by the
[AMR package](https://amr-for-r.org).

**Note on the data.** I don't have real hospital data, so this runs on `example_isolates`,
a dataset that ships with the AMR package. It holds 2,000 made-up isolates from 2002 to 2017.
To give the cleaning step something to do, `R/00_fetch.R` messes it up on purpose: organism
names spelled several ways, results written as "s", "Resistant" or "SENSITIVE", and some dates
in dd/mm/yyyy. None of the numbers describe a real hospital or country.

## Results

- [`reports/summary.md`](reports/summary.md): the main numbers and what I read from them
- [`reports/cleaning_log.md`](reports/cleaning_log.md): what cleaning changed, and which
  organism names a person should double-check
- `reports/tables/`: antibiograms, resistance by period and trend tests as CSV

![Resistance by period](reports/figures/resistance_trend.png)
![Antibiogram](reports/figures/antibiogram.png)

## How to run it

You need R (tested on 4.6.1). From the project folder:

```sh
Rscript -e 'renv::restore()'
Rscript R/run_all.R
```

This cleans the data, runs the analysis, writes the reports and then runs the tests.
[`docs/WORKFLOW.md`](docs/WORKFLOW.md) walks through each step, why it's done that way, and
how to point the scripts at a different export.

## What's in here

```
R/00_fetch.R      exports the example data and adds the deliberate mess
R/01_clean.R      organism names, S/I/R results, dates, duplicates
R/02_analyse.R    first isolates, antibiograms, resistance by period, trend test, charts
R/03_report.R     writes reports/summary.md and checks its claims against the numbers
R/run_all.R       runs all of the above, then the tests
data/raw/                  the messy input
data/organism_lookup.csv   manual fixes for names the AMR package matches wrongly
data/processed/            the cleaned isolates
tests/testthat/            tests for the cleaning step
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
- I used `example_isolates` rather than the package's `WHONET` example because `WHONET` only
  has 500 isolates spread over two countries and 16 years, which is too few for trends.

## Licence and credits

The code is released under the GNU General Public License v2 (see [`LICENSE`](LICENSE)).

The data in `data/` is derived from the `example_isolates` dataset in the AMR package, which
is also GPL-2. Intrinsic resistance rules come from EUCAST Expected Resistant Phenotypes
v1.2 (2023), as bundled in the AMR package.

If you use the AMR package, please cite it:

> Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C (2022). AMR: An R Package
> for Working with Antimicrobial Resistance Data. *Journal of Statistical Software*, 104(3),
> 1-31. [doi:10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)
