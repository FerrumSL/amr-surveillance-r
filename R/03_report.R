# Write <cfg$reports>/summary.md from `results` (made by R/02_analyse.R).
#
# Every number comes from the analysis. For the demo data the interpretation sentences are
# written by hand, so each claim they make is checked with stopifnot() first: if the data no
# longer support a claim, the run stops here and the text must be rewritten. For any other
# data the summary holds the numbers only; interpretation is left to a person.

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
})
if (!exists("results")) source("R/02_analyse.R")

p <- results$params
bp <- results$by_period
tr <- results$trend

pct <- function(x) sprintf("%.0f%%", x)
ci_txt <- function(lo, hi) sprintf("95%% CI %.0f-%.0f%%", lo, hi)
row_txt <- function(r) sprintf("%s (%s, n = %d) in %s", pct(r$pct_resistant),
                               ci_txt(r$ci_low, r$ci_high), r$n_tested, r$period)
pair_line <- function(pair_name) {
  d <- filter(bp, pair == pair_name, enough_isolates) |> arrange(period)
  t <- filter(tr, pair == pair_name)
  if (nrow(d) == 0) {
    return(sprintf("- **%s**: too few tested isolates per period for a percentage (%d first isolates).",
                   pair_name, t$n_isolates))
  }
  if (nrow(d) == 1) {
    return(sprintf("- **%s**: %s, the only period with at least %d tested isolates; %d first isolates.",
                   pair_name, row_txt(d), p$minimum, t$n_isolates))
  }
  # No trend is fitted when resistance never varies or the data cover a single year.
  if (is.null(t$p_value) || is.na(t$p_value)) {
    return(sprintf("- **%s**: %s, to %s. No trend test (resistance did not vary); %d first isolates.",
                   pair_name, row_txt(d[1, ]), row_txt(d[nrow(d), ]), t$n_isolates))
  }
  sprintf("- **%s**: %s, to %s. Odds ratio per year %.2f (95%% CI %.2f-%.2f), p = %s; %d first isolates.",
          pair_name, row_txt(d[1, ]), row_txt(d[nrow(d), ]),
          t$odds_ratio_per_year, t$or_ci_low, t$or_ci_high, format.pval(t$p_value, digits = 2),
          t$n_isolates)
}

alerts <- read_csv(file.path(p$reports, "tables/alerts.csv"), col_types = cols(.default = col_character()))
alert_lines <- if (nrow(alerts) == 0) "- No alerts." else
  count(alerts, alert) |> with(sprintf("- %s: %d isolates (%d patients)", alert, n,
                                        tapply(alerts$patient, alerts$alert, n_distinct)[alert]))

interpretation <- c("## Interpretation", "", "*Not generated for this data. Write it by hand after reading the tables.*")
if (isTRUE(p$demo)) {
  get_row <- function(ab_code, which_period) {
    filter(bp, ab == ab_code, period == which_period)
  }
  pct_by_period <- function(ab_code) filter(bp, ab == ab_code) |> arrange(period) |> pull(pct_resistant)
  first_period <- sort(unique(bp$period))[1]
  last_period <- sort(unique(bp$period), decreasing = TRUE)[1]
  amc <- filter(tr, ab == "AMC")
  cip <- filter(tr, ab == "CIP")
  sxt <- filter(tr, ab == "SXT")

  # The text below is written for exactly these three E. coli pairs.
  stopifnot(
    identical(sort(tr$ab), c("AMC", "CIP", "SXT")), all(tr$mo == "B_ESCHR_COLI"),
    all(bp$enough_isolates)
  )
  # Claims made in the interpretation below.
  stopifnot(
    amc$p_value < 0.05, amc$or_ci_low > 1,
    get_row("AMC", last_period)$ci_low > max(filter(bp, ab == "AMC", period != last_period)$pct_resistant),
    cip$p_value < 0.05, cip$or_ci_low > 1,
    get_row("CIP", first_period)$ci_high > get_row("CIP", last_period)$ci_low,
    cip$p_value > 0.05 / nrow(tr),
    sxt$p_value > 0.05,
    all(pct_by_period("SXT") > pmax(pct_by_period("AMC"), pct_by_period("CIP"))),
    # CoNS reported as a group until 2012, then almost only as species.
    with(results$first_by_year, sum(n_cons[year <= 2012]) > 150 && sum(n_cons[year >= 2014]) < 10)
  )
  amc_early <- filter(bp, ab == "AMC", period != last_period)$pct_resistant

  interpretation <- c(
    "## Interpretation",
    "",
    "*Hand-written; each claim is checked against the numbers by `R/03_report.R`.*",
    "",
    sprintf(paste(
      "E. coli resistance to amoxicillin/clavulanic acid stayed between %s and %s in the first %d",
      "periods and rose to %s in %s; the logistic trend test supports an increase (odds ratio %.2f per year).",
      "A straight-line trend on the log-odds scale is a simplification here, because most of the rise",
      "is in the last period."),
      pct(min(amc_early)), pct(max(amc_early)), length(amc_early),
      pct(get_row("AMC", last_period)$pct_resistant), last_period, amc$odds_ratio_per_year),
    "",
    sprintf(paste(
      "Ciprofloxacin resistance rose from %s in %s to %s in %s (odds ratio %.2f per year, p = %s).",
      "The evidence is weak: the confidence intervals of the first and last periods overlap, and the",
      "p-value would not pass a Bonferroni correction for the %d pairs tested (threshold %.3f)."),
      pct(get_row("CIP", first_period)$pct_resistant), first_period,
      pct(get_row("CIP", last_period)$pct_resistant), last_period, cip$odds_ratio_per_year,
      format.pval(cip$p_value, digits = 2), nrow(tr), 0.05 / nrow(tr)),
    "",
    sprintf(paste(
      "Trimethoprim/sulfamethoxazole resistance was the highest of the three pairs (%s-%s by",
      "period) and showed no clear trend (p = %s)."),
      pct(min(pct_by_period("SXT"))), pct(max(pct_by_period("SXT"))),
      format.pval(sxt$p_value, digits = 2)),
    "",
    paste(
      sprintf("Limits: each period holds only %d-%d tested E. coli isolates, three pairs were tested,",
              min(bp$n_tested), max(bp$n_tested)),
      "and the data come from one simulated setting. The trend test treats isolates from the same",
      "patient in different years as independent, which makes its p-values somewhat too small. In the same",
      "data, coagulase-negative staphylococci stop being reported as \"CoNS\" in 2013 and appear",
      "as species (S. epidermidis, S. hominis) instead, a change in laboratory practice that would",
      "distort any trend by species. See `reports/figures/isolates_per_year.png`.")
  )
}

lines <- c(
  "# Summary of resistance trends",
  "",
  sprintf("Generated by `R/03_report.R` for **%s**.", p$name),
  if (isTRUE(p$demo)) c("**Practice analysis on the AMR package's bundled `example_isolates` data",
                        "(randomised, fictitious isolates). Not real surveillance data.**"),
  "",
  "## Key figures",
  "",
  sprintf("- %d isolates after cleaning; %d first isolates (%s, %d-day episodes) from %d-%d.",
          results$n_isolates, results$n_first, p$first_isolate_method, p$episode_days,
          min(p$years), max(p$years)),
  sprintf("- Raw MIC and disk values interpreted with %s; I counted as %s. Percentages are only shown with at least %d tested isolates.",
          p$guideline, if (p$guideline_family == "CLSI") "not susceptible" else "susceptible", p$minimum),
  vapply(unique(bp$pair), pair_line, character(1)),
  "",
  "## Infection-control alerts",
  "",
  "All isolates, not only first isolates. Intrinsic resistance is excluded.",
  "",
  alert_lines,
  "",
  interpretation,
  "",
  "## Outputs",
  "",
  sprintf("- Printable antibiogram: `%s`", file.path(p$reports, "antibiogram.pdf")),
  sprintf("- Antibiogram tables: `%s/tables/antibiogram_gram_negative.csv`, `antibiogram_gram_positive.csv`", p$reports),
  sprintf("- Resistance by period: `%s/tables/resistance_by_period.csv`, chart `figures/resistance_trend.png`", p$reports),
  sprintf("- Trend tests: `%s/tables/trend_tests.csv`", p$reports),
  sprintf("- Alert list: `%s/tables/alerts.csv`", p$reports),
  sprintf("- Cleaning log: `%s/cleaning_log.md`", p$reports)
)
writeLines(lines, file.path(p$reports, "summary.md"))
message("Wrote ", file.path(p$reports, "summary.md"))
