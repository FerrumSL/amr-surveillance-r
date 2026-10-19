# Antibiograms, resistance by period and a simple trend test on first isolates.
#
# Input:  data/processed/isolates_clean.csv
# Output: reports/tables/*.csv, reports/figures/*.png, and `results` (used by R/03_report.R)

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ---- Parameters: change these for a new dataset --------------------------------------------
params <- list(
  guideline = "EUCAST",           # EUCAST: I counts as susceptible. CLSI: I counts as resistant.
  first_isolate_method = "episode-based",
  episode_days = 365,             # first isolate per patient and species per 365 days (CLSI M39)
  years = 2002:2017,              # analysis window
  period_length = 4,              # years per period in the trend table and chart
  minimum = 30,                   # no percentage is reported for fewer tested isolates
  # Gram-negative and Gram-positive panels are kept apart so that drugs an organism is
  # intrinsically resistant to (e.g. vancomycin for E. coli) do not show up as "resistance".
  antibiogram_gram_negative = c("AMX", "AMC", "CXM", "CTX", "CAZ", "MEM", "GEN", "CIP", "SXT", "NIT"),
  antibiogram_gram_positive = c("PEN", "FLC", "ERY", "CLI", "GEN", "CIP", "SXT", "TCY", "VAN", "RIF"),
  pairs = tribble(
    ~mo,            ~ab,
    "B_ESCHR_COLI", "CIP",
    "B_ESCHR_COLI", "AMC",
    "B_ESCHR_COLI", "SXT"
  )
)
options(AMR_guideline = params$guideline)

# ---- Data ------------------------------------------------------------------------------------
isolates <- read_csv("data/processed/isolates_clean.csv",
                     col_types = cols(.default = col_character()), na = "")
ab_cols <- setdiff(names(isolates), c("date", "patient", "age", "gender", "ward", "mo"))
isolates <- isolates |>
  mutate(date = as.Date(date), age = as.integer(age), mo = as.mo(mo, info = FALSE),
         across(all_of(ab_cols), as.sir))

# First isolates are picked from all dates, then cut to the analysis window.
first <- isolates |>
  filter_first_isolate(col_date = "date", col_patient_id = "patient", col_mo = "mo",
                       method = params$first_isolate_method,
                       episode_days = params$episode_days, info = FALSE) |>
  mutate(year = as.integer(format(date, "%Y"))) |>
  filter(year %in% params$years) |>
  mutate(period_start = min(params$years) +
           (year - min(params$years)) %/% params$period_length * params$period_length,
         period = paste0(period_start, "-",
                         pmin(period_start + params$period_length - 1, max(params$years))))
isolates <- filter(isolates, as.integer(format(date, "%Y")) %in% params$years)

# ---- Antibiograms ----------------------------------------------------------------------------
# Organisms with at least `minimum` first isolates; cells with fewer tested isolates stay empty.
# Intrinsic resistance (EUCAST Expected Resistant Phenotypes, e.g. K. pneumoniae to
# amoxicillin) is blanked: it is expected biology, not a surveillance finding.
make_antibiogram <- function(data, antibiotics) {
  common <- data |> count(mo) |> filter(n >= params$minimum) |> pull(mo)
  data <- filter(data, mo %in% common)
  for (ab in antibiotics) {
    data[[ab]][mo_is_intrinsic_resistant(data$mo, ab)] <- NA
  }
  antibiogram(data, antimicrobials = antibiotics,
              minimum = params$minimum, guideline = params$guideline, info = FALSE)
}
abg_neg <- make_antibiogram(filter(first, mo_is_gram_negative(mo)), params$antibiogram_gram_negative)
abg_pos <- make_antibiogram(filter(first, mo_is_gram_positive(mo)), params$antibiogram_gram_positive)

# ---- Resistance by period, with exact (Clopper-Pearson) 95% CIs --------------------------------
# I count as resistant only under CLSI, matching AMR::resistance(). The stopifnot() below
# checks that the hand-built flag gives the same proportion as the package.
is_resistant <- function(x) {
  if (params$guideline == "CLSI") x %in% c("R", "I") else x == "R"
}

pair_data <- bind_rows(Map(function(mo_code, ab_code) {
  first |>
    filter(mo == mo_code, !is.na(.data[[ab_code]])) |>
    transmute(mo = mo_code, ab = ab_code, year, period, result = .data[[ab_code]])
}, params$pairs$mo, params$pairs$ab)) |>
  mutate(resistant = is_resistant(result),
         organism = mo_shortname(mo, language = NULL),
         antibiotic = ab_name(ab, language = NULL),
         pair = paste(organism, "vs", tolower(antibiotic)))

by_period <- pair_data |>
  group_by(pair, organism, antibiotic, mo, ab, period) |>
  summarise(n_tested = n(), n_resistant = sum(resistant),
            check = resistance(result, minimum = 0, guideline = params$guideline),
            .groups = "drop") |>
  mutate(pct_resistant = 100 * n_resistant / n_tested)
stopifnot(isTRUE(all.equal(by_period$check, by_period$n_resistant / by_period$n_tested)))

ci <- mapply(function(k, n) binom.test(k, n)$conf.int, by_period$n_resistant, by_period$n_tested)
by_period <- by_period |>
  mutate(ci_low = 100 * ci[1, ], ci_high = 100 * ci[2, ],
         enough_isolates = n_tested >= params$minimum,
         across(c(pct_resistant, ci_low, ci_high), \(v) ifelse(enough_isolates, v, NA_real_))) |>
  select(-check)

# ---- Trend test: logistic regression of resistance on calendar year --------------------------
# One row per first isolate. Odds ratio per year > 1 means resistance odds rise over time.
# Profile-likelihood CI and likelihood-ratio p-value.
trend <- pair_data |>
  group_by(pair, organism, antibiotic, mo, ab) |>
  group_modify(function(d, key) {
    fit <- glm(resistant ~ year, family = binomial, data = d)
    ci <- suppressMessages(confint(fit, "year"))
    tibble(n_isolates = nrow(d), n_resistant = sum(d$resistant),
           first_year = min(d$year), last_year = max(d$year),
           odds_ratio_per_year = exp(coef(fit)[["year"]]),
           or_ci_low = exp(ci[[1]]), or_ci_high = exp(ci[[2]]),
           p_value = anova(fit, test = "LRT")["year", "Pr(>Chi)"])
  }) |>
  ungroup()

# ---- Tables ----------------------------------------------------------------------------------
dir.create("reports/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("reports/figures", showWarnings = FALSE, recursive = TRUE)
write_csv(as.data.frame(abg_neg), "reports/tables/antibiogram_gram_negative.csv", na = "")
write_csv(as.data.frame(abg_pos), "reports/tables/antibiogram_gram_positive.csv", na = "")
write_csv(attr(abg_neg, "long_numeric"), "reports/tables/antibiogram_gram_negative_long.csv", na = "")
write_csv(attr(abg_pos, "long_numeric"), "reports/tables/antibiogram_gram_positive_long.csv", na = "")
write_csv(select(by_period, -mo, -ab), "reports/tables/resistance_by_period.csv", na = "")
write_csv(select(trend, -mo, -ab), "reports/tables/trend_tests.csv", na = "")

# ---- Figures ---------------------------------------------------------------------------------
fig_trend <- by_period |>
  filter(enough_isolates) |>
  ggplot(aes(period, pct_resistant, group = pair)) +
  geom_line(colour = "grey60") +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), colour = "#b2182b") +
  geom_text(aes(y = ci_high, label = paste0("n=", n_tested)), vjust = -0.6, size = 3) +
  facet_wrap(~pair, ncol = 1) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Resistance by period (first isolates)",
       subtitle = paste0("Points: % resistant (", params$guideline,
                         "). Bars: exact 95% CI. n = isolates tested."),
       x = NULL, y = "% resistant") +
  theme_minimal(base_size = 11)
ggsave("reports/figures/resistance_trend.png", fig_trend, width = 7, height = 7, dpi = 150)

heat <- bind_rows(
  mutate(attr(abg_neg, "long_numeric"), panel = "Gram-negative"),
  mutate(attr(abg_pos, "long_numeric"), panel = "Gram-positive")
) |>
  filter(!is.na(coverage))
fig_abg <- ggplot(heat, aes(ab, mo, fill = 100 * coverage)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = paste0(round(100 * coverage), "%\nn=", n_tested)), size = 2.6) +
  scale_fill_gradient(low = "#d7301f", high = "#fef0d9", limits = c(0, 100)) +
  facet_wrap(~panel, ncol = 1, scales = "free") +
  labs(title = "Antibiogram: % susceptible (first isolates)",
       subtitle = paste0("Only cells with at least ", params$minimum, " tested isolates. ",
                         params$guideline, " interpretation."),
       x = NULL, y = NULL, fill = "% S") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
ggsave("reports/figures/antibiogram.png", fig_abg, width = 8, height = 6.5, dpi = 150)

fig_n <- first |>
  mutate(organism = mo_shortname(mo, language = NULL),
         organism = ifelse(organism %in% names(sort(table(organism), decreasing = TRUE))[1:5],
                           organism, "Other")) |>
  count(year, organism) |>
  ggplot(aes(year, n, fill = organism)) +
  geom_col() +
  labs(title = "First isolates per year", subtitle = "Five most common organisms",
       x = NULL, y = "First isolates", fill = NULL) +
  theme_minimal(base_size = 11)
ggsave("reports/figures/isolates_per_year.png", fig_n, width = 7, height = 4, dpi = 150)

first_by_year <- first |>
  group_by(year) |>
  summarise(n = n(), n_cons = sum(mo == "B_STPHY_CONS"))

results <- list(params = params, n_isolates = nrow(isolates), n_first = nrow(first),
                first_by_year = first_by_year,
                abg_neg = abg_neg, abg_pos = abg_pos, by_period = by_period, trend = trend)
message("Wrote reports/tables/ and reports/figures/")
