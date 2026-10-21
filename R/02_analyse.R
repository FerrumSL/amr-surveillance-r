# Antibiograms, resistance by period and a simple trend test on first isolates.
#
# Input:  cfg$processed (settings in config/*.R)
# Output: <cfg$reports>/tables/*.csv, figures/*.png, antibiogram.pdf, and `results` (used by R/03_report.R)

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

if (!exists("cfg")) source("config/demo.R")
params <- cfg
# EUCAST counts I ("susceptible, increased exposure") as susceptible. Under CLSI, I is counted
# as not susceptible here, as the AMR package does.
params$guideline_family <- sub(" .*", "", cfg$guideline)
report_path <- function(...) file.path(params$reports, ...)
options(AMR_guideline = params$guideline_family)

# ---- Data ------------------------------------------------------------------------------------
isolates <- read_csv(params$processed, col_types = cols(.default = col_character()), na = "")
ab_cols <- setdiff(names(isolates), c("date", "patient", "age", "gender", "ward", "specimen", "mo"))
isolates <- isolates |>
  mutate(date = as.Date(date), age = as.integer(age), mo = as.mo(mo, info = FALSE),
         across(all_of(ab_cols), as.sir))
data_years <- as.integer(format(isolates$date, "%Y"))
if (is.null(params$years)) params$years <- min(data_years):max(data_years)
if (is.null(params$antibiogram_years)) params$antibiogram_years <- params$years

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
isolates <- filter(isolates, data_years %in% params$years)

# ---- Antibiograms ----------------------------------------------------------------------------
# Organisms with at least `minimum` first isolates; cells with fewer tested isolates stay empty.
# Intrinsic resistance (EUCAST Expected Resistant Phenotypes, e.g. K. pneumoniae to
# amoxicillin) is blanked: it is expected biology, not a surveillance finding.
make_antibiogram <- function(data, antibiotics) {
  antibiotics <- intersect(antibiotics, names(data))
  common <- data |> count(mo) |> filter(n >= params$minimum) |> pull(mo)
  data <- filter(data, mo %in% common)
  for (ab in antibiotics) {
    data[[ab]][mo_is_intrinsic_resistant(data$mo, ab)] <- NA
  }
  antibiogram(data, antimicrobials = antibiotics,
              minimum = params$minimum, guideline = params$guideline_family, info = FALSE)
}
abg_first <- filter(first, year %in% params$antibiogram_years)
abg_neg <- make_antibiogram(filter(abg_first, mo_is_gram_negative(mo)), params$antibiogram_gram_negative)
abg_pos <- make_antibiogram(filter(abg_first, mo_is_gram_positive(mo)), params$antibiogram_gram_positive)

# ---- Resistance by period, with exact (Clopper-Pearson) 95% CIs --------------------------------
# I count as resistant only under CLSI, matching AMR::resistance(). The stopifnot() below
# checks that the hand-built flag gives the same proportion as the package.
is_resistant <- function(x) {
  if (params$guideline_family == "CLSI") x %in% c("R", "I") else x == "R"
}

# A species code also matches its subspecies (B_KLBSL_PNMN covers B_KLBSL_PNMN_PNMN).
pair_data <- bind_rows(Map(function(mo_code, ab_code) {
  if (!ab_code %in% names(first)) return(NULL)
  first |>
    # NI (not interpretable) is a result but not a test outcome, so it stays out of n_tested.
    filter(mo == mo_code | startsWith(as.character(mo), paste0(mo_code, "_")),
           .data[[ab_code]] %in% c("S", "SDD", "I", "R")) |>
    transmute(mo = mo_code, ab = ab_code, year, period, result = .data[[ab_code]])
}, params$pairs$mo, params$pairs$ab)) |>
  mutate(resistant = is_resistant(result),
         organism = mo_shortname(mo, language = NULL),
         antibiotic = ab_name(ab, language = NULL),
         pair = paste(organism, "vs", tolower(antibiotic)))

by_period <- pair_data |>
  group_by(pair, organism, antibiotic, mo, ab, period) |>
  summarise(n_tested = n(), n_resistant = sum(resistant),
            check = resistance(result, minimum = 0, guideline = params$guideline_family),
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
    out <- tibble(n_isolates = nrow(d), n_resistant = sum(d$resistant),
                  first_year = min(d$year), last_year = max(d$year))
    # No trend can be fitted to a single year or to all-resistant / all-susceptible data.
    if (n_distinct(d$year) < 2 || n_distinct(d$resistant) < 2) return(out)
    fit <- glm(resistant ~ year, family = binomial, data = d)
    ci <- suppressMessages(confint(fit, "year"))
    mutate(out, odds_ratio_per_year = exp(coef(fit)[["year"]]),
           or_ci_low = exp(ci[[1]]), or_ci_high = exp(ci[[2]]),
           p_value = anova(fit, test = "LRT")["year", "Pr(>Chi)"])
  }) |>
  ungroup()

# ---- Tables ----------------------------------------------------------------------------------
dir.create(report_path("tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(report_path("figures"), showWarnings = FALSE, recursive = TRUE)
write_csv(as.data.frame(abg_neg), report_path("tables/antibiogram_gram_negative.csv"), na = "")
write_csv(as.data.frame(abg_pos), report_path("tables/antibiogram_gram_positive.csv"), na = "")
write_csv(attr(abg_neg, "long_numeric"), report_path("tables/antibiogram_gram_negative_long.csv"), na = "")
write_csv(attr(abg_pos, "long_numeric"), report_path("tables/antibiogram_gram_positive_long.csv"), na = "")
write_csv(select(by_period, -mo, -ab), report_path("tables/resistance_by_period.csv"), na = "")
write_csv(select(trend, -mo, -ab), report_path("tables/trend_tests.csv"), na = "")

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
ggsave(report_path("figures/resistance_trend.png"), fig_trend, width = 7, height = 7, dpi = 150)

# One page for clinicians, laid out as CLSI M39 asks: organisms with their number of first
# isolates, % susceptible per drug, blank where fewer than `minimum` were tested.
heat <- bind_rows(
  mutate(attr(abg_neg, "long_numeric"), panel = "Gram-negative"),
  mutate(attr(abg_pos, "long_numeric"), panel = "Gram-positive")
) |>
  filter(!is.na(coverage)) |>
  mutate(organism = paste0(mo, " (", n_total, ")"),
         pct_s = round(100 * coverage),  # the colour band uses the printed number
         band = cut(pct_s, c(-Inf, 60, 90, Inf), right = FALSE,
                    labels = c("< 60%", "60-89%", ">= 90%")))
# Some drugs have no S category in the guideline: colistin under CLSI since 2020 is I or R only,
# stored by the AMR package as an S breakpoint of 0.0001. "0% susceptible" would read as total
# resistance, so those cells say n/a. Drugs that do have an S breakpoint keep their 0%.
no_s <- clinical_breakpoints |>
  filter(guideline == params$guideline, type == "human") |>
  mutate(ab = as.character(ab)) |>
  group_by(ab) |>
  summarise(no_s = all(breakpoint_S < 0.001)) |>
  filter(no_s) |>
  pull(ab)
heat <- mutate(heat, no_s = as.character(as.ab(ab, info = FALSE)) %in% no_s & pct_s == 0,
               label = ifelse(no_s, "n/a*", pct_s),
               band = replace(band, no_s, NA))
abg_span <- range(params$antibiogram_years)
fig_abg <- ggplot(heat, aes(ab, organism, fill = band)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.4) +
  scale_x_discrete(labels = \(x) sub("/", "/\n", x)) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("< 60%" = "#f4a582", "60-89%" = "#fde9b5", ">= 90%" = "#d1e5f0"),
                    drop = FALSE, na.value = "grey92", na.translate = FALSE) +
  facet_wrap(~panel, ncol = 1, scales = "free") +
  labs(title = paste0("Cumulative antibiogram, ", abg_span[1],
                      if (abg_span[2] > abg_span[1]) paste0("-", abg_span[2])),
       subtitle = paste0(params$name, "\n% susceptible. Organism (number of first isolates)."),
       caption = paste0("First isolate per patient and species per ", params$episode_days,
                        " days (CLSI M39). Raw MIC and disk values: ", params$guideline,
                        " breakpoints; S/I/R reported by the lab kept as reported.\n",
                        "Blank: fewer than ", params$minimum, " isolates tested, intrinsic resistance, or not tested. ",
                        "n/a*: the guideline has no susceptible category for this drug (see the CSV tables for % intermediate).\n",
                        "Generated ", format(Sys.Date()), " by R/02_analyse.R."),
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom",
        plot.caption = element_text(hjust = 0))
ggsave(report_path("figures/antibiogram.png"), fig_abg, width = 11.69, height = 8.27, dpi = 150)
ggsave(report_path("antibiogram.pdf"), fig_abg, width = 11.69, height = 8.27)  # A4 landscape

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
ggsave(report_path("figures/isolates_per_year.png"), fig_n, width = 7, height = 4, dpi = 150)

first_by_year <- first |>
  group_by(year) |>
  summarise(n = n(), n_cons = sum(mo == "B_STPHY_CONS"))

results <- list(params = params, n_isolates = nrow(isolates), n_first = nrow(first),
                first_by_year = first_by_year,
                abg_neg = abg_neg, abg_pos = abg_pos, by_period = by_period, trend = trend)
message("Wrote ", params$reports, "/tables, figures and antibiogram.pdf")
