# Antibiograms, resistance by period and a simple trend test on first isolates.
#
# Input:  cfg$processed (settings in config/*.R)
# Output: <cfg$reports>/tables/*.csv, figures/*.png, antibiogram.pdf, and `results` (used by R/03_report.R)

suppressPackageStartupMessages({
  library(AMR)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(grid)
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
ab_cols <- setdiff(names(isolates), c("date", "patient", "age", "gender", "ward", "specimen", "dtest", "mo"))
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
# Palette, theme and the print checks live in R/plot_style.R. The checks run on every build, so
# a palette edit that breaks black-and-white print or text contrast stops the run here.
source("R/plot_style.R")
invisible(check_palette(band_fill))
category_all <- c(category_fill, Other = other_fill)
invisible(check_palette(category_all[order(-lstar(category_all))], text = FALSE))
gray_dir <- file.path("reports", "_grayscale_check")
dir.create(gray_dir, showWarnings = FALSE, recursive = TRUE)
run_name <- basename(params$reports)

pairs_shown <- n_distinct(filter(by_period, enough_isolates)$pair)
fig_trend <- by_period |>
  filter(enough_isolates) |>
  ggplot(aes(period, pct_resistant, group = pair)) +
  geom_line(colour = band_fill[["60-79%"]], linewidth = 0.7) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), colour = band_fill[["20-39%"]], linewidth = 0.8) +
  geom_point(colour = band_fill[["20-39%"]], size = 2.6) +
  geom_text(aes(y = ci_high, label = paste0("n = ", n_tested)), vjust = -0.9, size = 2.7, colour = ink_2) +
  facet_wrap(~pair, ncol = 1) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), labels = \(x) paste0(x, "%"),
                     expand = c(0, 0)) +
  labs(title = paste0("Resistance by ", params$period_length, "-year period"),
       subtitle = paste0(params$name, "\nFirst isolates. Point: % resistant. Bar: exact 95% CI."),
       caption = paste0(params$guideline, " interpretation. Periods with fewer than ", params$minimum,
                        " tested isolates are not shown."),
       x = NULL, y = "% resistant") +
  theme_lab(10) +
  theme(panel.spacing.y = unit(14, "pt"),
        axis.title.y = element_text(face = "bold", size = 11, colour = ink, margin = margin(r = 8)))
ggsave(report_path("figures/resistance_trend.png"), fig_trend, width = 7,
       height = 1.6 + 2 * max(pairs_shown, 1), dpi = 150, bg = "white")

# ---- Printable antibiogram -------------------------------------------------------------------
# One A4 page for clinicians, laid out as CLSI M39 asks: organisms with their number of first
# isolates, % susceptible per drug, drugs grouped by class. A cell is one of:
#   value  % susceptible with n tested (italic with a dagger for 30-49 tested)
#   n/a*   the guideline has no S category for the drug (colistin under CLSI)
#   IR     intrinsic resistance, so not reported
#   -      tested, but fewer than `minimum` isolates
#   blank  not tested
drug_class <- function(code) {
  g <- ab_group(code, language = NULL)
  case_when(grepl("penicillin", g, ignore.case = TRUE) ~ "Penicillins",
            grepl("^Cephalosporin", g) ~ "Cephalosporins",
            g %in% c("Carbapenems", "Aminoglycosides", "Fluoroquinolones", "Glycopeptides") ~ g,
            TRUE ~ "Other")
}
class_order <- c("Penicillins", "Cephalosporins", "Carbapenems", "Aminoglycosides",
                 "Fluoroquinolones", "Glycopeptides", "Other")
class_short <- c(Penicillins = "Peni-\ncillins", Cephalosporins = "Cephalo-\nsporins",
                 Carbapenems = "Carba-\npenems", Aminoglycosides = "Amino-\nglycosides",
                 Fluoroquinolones = "Fluoro-\nquinolones", Glycopeptides = "Glyco-\npeptides",
                 Other = "Other")

# Drugs with no S category in the guideline: stored by the AMR package as an S breakpoint of
# 0.0001. "0% susceptible" would read as total resistance, so those cells say n/a.
no_s <- clinical_breakpoints |>
  filter(guideline == params$guideline, type == "human") |>
  mutate(ab = as.character(ab)) |>
  group_by(ab) |>
  summarise(no_s = all(breakpoint_S < 0.001)) |>
  filter(no_s) |>
  pull(ab)

# Every organism x drug cell of one panel, with its state. The values come unchanged from the
# antibiogram() table; the grid only adds what that table leaves out (IR, too few, not tested).
# Staphylococcus species with erythromycin-R, clindamycin-S isolates and no D-test result.
dtest_missing <- if (all(c("ERY", "CLI") %in% names(abg_first))) {
  unique(as.character(abg_first$mo[abg_first$ERY %in% "R" & abg_first$CLI %in% "S" & is.na(abg_first$dtest) &
                                     mo_genus(abg_first$mo, language = NULL) == "Staphylococcus"]))
} else character()
abg_cells <- function(data, antibiotics, abg) {
  antibiotics <- intersect(antibiotics, names(data))
  common <- data |> count(mo) |> filter(n >= params$minimum)
  if (nrow(common) == 0 || length(antibiotics) == 0) return(NULL)
  data <- filter(data, mo %in% common$mo)
  long <- attr(abg, "long_numeric")
  cells <- bind_rows(lapply(antibiotics, function(ab) {
    data |>
      group_by(mo) |>
      summarise(tested = sum(!is.na(.data[[ab]])), .groups = "drop") |>
      mutate(code = ab, intrinsic = mo_is_intrinsic_resistant(mo, ab))
  })) |>
    left_join(common, by = "mo") |>
    mutate(organism = mo_shortname(mo, language = NULL), ab = ab_name(code, language = NULL)) |>
    left_join(select(long, organism = mo, ab, coverage, n_tested), by = c("organism", "ab")) |>
    mutate(pct_s = round(100 * coverage),
           state = case_when(!is.na(coverage) & code %in% no_s & pct_s == 0 ~ "no_s",
                             !is.na(coverage) ~ "value",
                             intrinsic ~ "IR",
                             tested > 0 ~ "few",
                             TRUE ~ "blank"),
           band = band_of(ifelse(state == "value", pct_s, NA)),
           wide = state == "value" & n_tested < 50,
           dtest_flag = state == "value" & code == "CLI" & mo_genus(mo, language = NULL) == "Staphylococcus" &
             mo %in% dtest_missing,
           class = factor(drug_class(code), class_order))
  # A drug with nothing to show for any organism (only dashes and blanks) is left off the page.
  shown <- unique(cells$code[cells$state %in% c("value", "no_s", "IR")])
  cells <- filter(cells, code %in% shown)
  # When every value in a row has 30-49 isolates, one dagger on the organism replaces them all.
  cells <- cells |>
    group_by(organism) |>
    mutate(row_dagger = any(state == "value") && all(wide[state == "value"]), wide = wide & !row_dagger) |>
    ungroup()
  # Columns ordered by class, then as listed in the config; a gap separates the classes.
  cols <- distinct(cells, class, code, ab) |>
    arrange(class, match(code, antibiotics)) |>
    mutate(x = row_number() + 0.35 * (as.integer(factor(class)) - 1))
  rows <- distinct(cells, organism, n) |> arrange(desc(organism)) |> mutate(y = row_number())
  cells |> left_join(select(cols, code, x), by = "code") |> left_join(select(rows, organism, y), by = "organism")
}
cells_neg <- abg_cells(filter(abg_first, mo_is_gram_negative(mo)), params$antibiogram_gram_negative, abg_neg)
cells_pos <- abg_cells(filter(abg_first, mo_is_gram_positive(mo)), params$antibiogram_gram_positive, abg_pos)

abg_panel <- function(cells, cell_w, cell_h, label_units, top_in, name_in, gray) {
  fill <- if (gray) setNames(to_gray(band_fill), names(band_fill)) else band_fill
  n_rows <- max(cells$y)
  top <- n_rows + 0.5 + top_in / cell_h
  vals <- filter(cells, state %in% c("value", "no_s")) |>
    mutate(fill = ifelse(state == "value", fill[as.character(band)], "#ffffff"),
           text = ifelse(state == "value", text_on(fill), ink),
           label = ifelse(state == "no_s", "n/a*", paste0(pct_s, ifelse(wide, "\u2020", ""),
                                                          ifelse(dtest_flag, "\u2021", ""))),
           face = ifelse(wide, "bold.italic", "bold"))
  cols <- distinct(cells, class, code, ab, x)
  brackets <- cols |>
    group_by(class) |>
    summarise(x0 = min(x) - 0.42, x1 = max(x) + 0.42, n = n(), .groups = "drop") |>
    mutate(label = ifelse(n == 1, class_short[as.character(class)], as.character(class)))
  bracket_y <- n_rows + 0.5 + (name_in + 0.16) / cell_h
  rows <- distinct(cells, organism, n, y, row_dagger)
  k <- min(1, cell_h / 0.45)  # text shrinks with the cells when many organisms share the page
  ggplot() +
    geom_tile(data = filter(vals, state == "value"), aes(x, y), fill = filter(vals, state == "value")$fill,
              width = 0.94, height = 0.9) +
    geom_tile(data = filter(vals, state == "no_s"), aes(x, y), fill = "#ffffff", colour = ink_2,
              linetype = "22", linewidth = 0.4, width = 0.86, height = 0.82) +
    geom_tile(data = filter(cells, state == "IR"), aes(x, y), fill = "#ffffff", colour = ink_3,
              linewidth = 0.3, width = 0.86, height = 0.82) +
    geom_text(data = vals, aes(x, y + 0.13, label = label), colour = vals$text, fontface = vals$face,
              size = max(3.9 * k, 3.4)) +  # >= 9.5 pt
    geom_text(data = filter(vals, state == "value"), aes(x, y - 0.26, label = paste("n =", n_tested)),
              colour = filter(vals, state == "value")$text, size = max(2.5 * k, 2.3)) +  # >= 6.5 pt
    geom_text(data = filter(vals, state == "no_s"), aes(x, y - 0.26, label = paste("n =", n_tested)),
              colour = ink_2, size = max(2.5 * k, 2.3)) +
    geom_text(data = filter(cells, state == "IR"), aes(x, y, label = "IR"), colour = ink_2, size = 3.2,
              fontface = "bold") +
    geom_text(data = filter(cells, state == "few"), aes(x, y, label = "—"), colour = ink_3, size = 4) +
    geom_text(data = cols, aes(x, n_rows + 0.5 + 0.08 / cell_h, label = sub("/", "/\n", ab)), angle = 90,
              hjust = 0, vjust = 0.5, size = 2.9, lineheight = 0.9, colour = ink) +
    geom_segment(data = brackets, aes(x = x0, xend = x1, y = bracket_y, yend = bracket_y),
                 colour = ink_2, linewidth = 0.35) +
    geom_segment(data = brackets, aes(x = x0, xend = x0, y = bracket_y, yend = bracket_y - 0.12 / cell_h),
                 colour = ink_2, linewidth = 0.35) +
    geom_segment(data = brackets, aes(x = x1, xend = x1, y = bracket_y, yend = bracket_y - 0.12 / cell_h),
                 colour = ink_2, linewidth = 0.35) +
    geom_text(data = brackets, aes((x0 + x1) / 2, bracket_y + 0.06 / cell_h, label = label),
              vjust = 0, size = 2.8, lineheight = 0.9, colour = ink_2) +
    geom_text(data = rows, aes(0.35, y, label = sprintf("italic('%s')~'(%d)%s'", organism, n,
                                                        ifelse(row_dagger, "\u2020", ""))),
              parse = TRUE, hjust = 1, size = 3.4, colour = ink) +
    coord_fixed(ratio = cell_h / cell_w, xlim = c(0.4 - label_units, max(cols$x) + 0.6),
                ylim = c(0.5, top), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

draw_antibiogram <- function(gray = FALSE) {
  grid.newpage()  # first: measuring text below needs an open page, and must not start a second
  page_w <- 11.69
  page_h <- 8.27
  head_h <- 1.15
  foot_h <- 1.45
  gap_h <- 0.42     # panel title and space between panels
  panels <- Filter(\(p) !is.null(p$cells), list(list(cells = cells_neg, title = "Gram-negative"),
                                               list(cells = cells_pos, title = "Gram-positive")))
  fill <- if (gray) setNames(to_gray(band_fill), names(band_fill)) else band_fill
  rows <- vapply(panels, \(p) max(p$cells$y), numeric(1))
  xmax <- vapply(panels, \(p) max(p$cells$x), numeric(1))
  label_in <- max(unlist(lapply(panels, \(p) vapply(unique(paste0(p$cells$organism, " (", p$cells$n, ")")),
    \(o) convertWidth(stringWidth(o), "in", valueOnly = TRUE), numeric(1))))) * 9.7 / 12 + 0.15
  # Height of the tallest vertical drug name (its longest line), then the bracket and class label.
  drug_lines <- unlist(strsplit(sub("/", "/\n", unique(unlist(lapply(panels, \(p) p$cells$ab)))), "\n"))
  name_in <- max(vapply(drug_lines, \(t) convertWidth(stringWidth(t), "in", valueOnly = TRUE), numeric(1))) * 8.2 / 12
  top_in <- name_in + 0.55
  cell_h <- min(0.58, (page_h - head_h - foot_h - length(panels) * (top_in + gap_h)) / max(sum(rows), 1))
  cell_h <- max(cell_h, 0.3)  # ponytail: more than ~12 organism rows would need a second page
  cell_w <- min(0.8, (page_w - 0.7 - label_in) / (max(xmax) + 0.6))
  label_units <- label_in / cell_w
  span <- range(params$antibiogram_years)
  title <- paste0("Cumulative antibiogram ", span[1], if (span[2] > span[1]) paste0("–", span[2]))

  grid.rect(gp = gpar(fill = "white", col = NA))
  heights <- c(head_h, rows * cell_h + top_in + gap_h, foot_h)
  pushViewport(viewport(layout = grid.layout(length(heights) + 1, 1,
    heights = unit(c(heights, 1), c(rep("in", length(heights)), "null")))))

  # Header: title, demo badge, subtitle with the breakpoint standard, and how to read a cell.
  pushViewport(viewport(layout.pos.row = 1))
  title_gp <- gpar(fontsize = 18, fontface = "bold", col = ink)
  grid.text(title, x = unit(0.3, "in"), y = unit(0.66, "npc"), just = "left", gp = title_gp)
  if (isTRUE(params$demo)) {
    bx <- 0.3 + convertWidth(grobWidth(textGrob(title, gp = title_gp)), "in", valueOnly = TRUE) + 0.2
    grid.roundrect(x = unit(bx, "in"), y = unit(0.66, "npc"), width = unit(1.95, "in"), height = unit(0.28, "in"),
                   just = "left", r = unit(0.05, "in"), gp = gpar(fill = ink, col = NA))
    grid.text("SYNTHETIC TEST DATA", x = unit(bx + 0.975, "in"), y = unit(0.66, "npc"),
              gp = gpar(fontsize = 9, fontface = "bold", col = "white"))
  }
  grid.text(paste0(params$name, "   ·   ", params$guideline, " breakpoints",
                   if (!is.null(params$antibiogram_subtitle)) paste0("   ·   ", params$antibiogram_subtitle)),
            x = unit(0.3, "in"), y = unit(0.33, "npc"), just = "left", gp = gpar(fontsize = 9.5, col = ink_2))
  grid.text("Large number: % susceptible\nn = isolates tested for that drug\n(n): first isolates of the organism",
            x = unit(1, "npc") - unit(0.3, "in"), y = unit(0.62, "npc"), just = "right",
            gp = gpar(fontsize = 8, col = ink_2, lineheight = 1.2))
  grid.lines(x = unit(c(0.3, page_w - 0.3), "in"), y = unit(0.1, "npc"), gp = gpar(col = rule, lwd = 1))
  popViewport()

  for (i in seq_along(panels)) {
    pushViewport(viewport(layout.pos.row = i + 1))
    grid.text(panels[[i]]$title, x = unit(0.3, "in"), y = unit(1, "npc") - unit(0.05, "in"),
              just = c("left", "top"), gp = gpar(fontsize = 11, fontface = "bold", col = ink))
    w <- (xmax[i] + 0.6 - 0.4) * cell_w + label_in
    h <- rows[i] * cell_h + top_in
    pushViewport(viewport(x = unit(0.3, "in"), y = unit(0.2, "in"), width = unit(w, "in"), height = unit(h, "in"),
                          just = c("left", "bottom")))
    print(abg_panel(panels[[i]]$cells, cell_w, cell_h, label_units, top_in, name_in, gray), newpage = FALSE)
    popViewport(2)
  }

  # Footer: the six bins, the three non-value cell types, and the notes.
  pushViewport(viewport(layout.pos.row = length(heights)))
  y_key <- unit(1, "npc") - unit(0.25, "in")
  grid.text("% susceptible", x = unit(0.3, "in"), y = y_key, just = "left",
            gp = gpar(fontsize = 8, fontface = "bold", col = ink))
  for (k in seq_along(fill)) {
    x <- 1.25 + (k - 1) * 0.85
    grid.rect(x = unit(x, "in"), y = y_key, width = unit(0.22, "in"), height = unit(0.15, "in"), just = "left",
              gp = gpar(fill = fill[[k]], col = NA))
    grid.text(names(fill)[k], x = unit(x + 0.28, "in"), y = y_key, just = "left", gp = gpar(fontsize = 8, col = ink))
  }
  key <- function(x, label, text, border, lty = "solid", face = "bold") {
    if (!is.na(border)) {
      grid.rect(x = unit(x, "in"), y = y_key, width = unit(0.32, "in"), height = unit(0.17, "in"), just = "left",
                gp = gpar(fill = "white", col = border, lty = lty, lwd = 0.7))
    }
    grid.text(label, x = unit(x + 0.16, "in"), y = y_key, gp = gpar(fontsize = 6.5, fontface = face, col = ink_2))
    grid.text(text, x = unit(x + 0.4, "in"), y = y_key, just = "left", gp = gpar(fontsize = 8, col = ink))
  }
  # Only the cell types and notes that occur on this page.
  all_cells <- bind_rows(lapply(panels, `[[`, "cells"))
  has <- function(st) any(all_cells$state == st)
  x <- 1.25 + 6 * 0.85 + 0.25
  if (has("no_s")) { key(x, "n/a*", "no S category", ink_2, lty = "22"); x <- x + 1.4 }
  if (has("IR")) { key(x, "IR", "intrinsic", ink_3); x <- x + 1.05 }
  if (has("few")) key(x, "\u2014", paste0("< ", params$minimum, " tested"), NA, face = "plain")
  no_s_drugs <- unique(all_cells$ab[all_cells$state == "no_s"])
  # Staphylococci that are erythromycin R and clindamycin S need a D-test (CLSI M100). Without
  # one, inducible clindamycin resistance cannot be excluded, so the page says how many.
  staph <- abg_first |> filter(mo_genus(mo, language = NULL) == "Staphylococcus")
  undetected <- if (all(c("ERY", "CLI") %in% names(staph)) && "CLI" %in% all_cells$code) {
    sum(staph$ERY %in% "R" & staph$CLI %in% "S" & is.na(staph$dtest))
  } else 0
  notes <- c(
    paste0("First isolate per patient and species per ", params$episode_days, " days (CLSI M39). ",
           "Raw MIC and disk values read with ", params$guideline,
           " breakpoints; S/I/R reported by the laboratory kept as reported."),
    paste(c(if (has("IR")) "IR: intrinsic resistance (EUCAST Expected Resistant Phenotypes v1.2), not reported.",
            if (has("few")) paste0("\u2014: fewer than ", params$minimum, " isolates tested."),
            if (has("blank")) "Blank: not tested.",
            if (any(all_cells$wide | all_cells$row_dagger))
              paste0("\u2020 ", params$minimum, "\u201349 isolates tested, wide uncertainty",
                     if (any(all_cells$row_dagger)) " (after an organism: every value in that row)", ".")),
          collapse = "   "),
    if (undetected > 0)
      paste0("\u2021 Clindamycin: ", undetected, " staphylococci are erythromycin R and clindamycin S with no ",
             "D-test result; inducible resistance cannot be excluded, so % susceptible may be too high."),
    if (length(no_s_drugs) > 0)
      paste0("n/a*: ", params$guideline, " has no susceptible category for ", paste(no_s_drugs, collapse = ", "),
             "; see the CSV tables for % intermediate."),
    paste0("Generated ", format(Sys.Date()), " by R/02_analyse.R. ",
           if (is.null(params$antibiogram_note)) paste("Use with local guidelines; an antibiogram guides empirical",
                                                       "therapy and does not replace the patient's own result.")
           else params$antibiogram_note))
  notes <- notes[nzchar(notes)]
  grid.text(paste(notes, collapse = "\n"), x = unit(0.3, "in"), y = y_key - unit(0.25, "in"), just = c("left", "top"),
            gp = gpar(fontsize = 8, col = ink_2, lineheight = 1.3))
  popViewport(2)
}
png(report_path("figures/antibiogram.png"), width = 11.69, height = 8.27, units = "in", res = 150, type = "cairo")
draw_antibiogram()
invisible(dev.off())
cairo_pdf(report_path("antibiogram.pdf"), width = 11.69, height = 8.27)  # A4 landscape; cairo for – † ·
draw_antibiogram()
invisible(dev.off())
png(file.path(gray_dir, paste0(run_name, "_antibiogram.png")), width = 11.69, height = 8.27, units = "in",
    res = 110, type = "cairo")
draw_antibiogram(gray = TRUE)
invisible(dev.off())

# ---- First isolates per year -----------------------------------------------------------------
# Four most common organisms and "Other"; hue and lightness both differ, so it survives grayscale.
by_org <- mutate(first, organism = mo_shortname(mo, language = NULL),
                 organism = ifelse(organism == "CoNS", "CoNS (not identified to species)", organism))
top4 <- by_org |> count(organism, sort = TRUE) |> slice_head(n = 4) |> pull(organism)
year_counts <- by_org |>
  mutate(organism = factor(ifelse(organism %in% top4, organism, "Other"), c(top4, "Other"))) |>
  count(year, organism)
n_years <- n_distinct(year_counts$year)
isolates_plot <- function(gray) {
  fills <- c(setNames(category_fill[seq_along(top4)], top4), Other = other_fill)
  if (gray) fills <- setNames(to_gray(fills), names(fills))
  ggplot(year_counts, aes(year, n, fill = organism)) +
    geom_col(colour = "white", linewidth = 0.3, width = 0.8) +
    scale_fill_manual(values = fills) +
    scale_x_continuous(breaks = sort(unique(year_counts$year)), expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "First isolates per year", subtitle = paste0(params$name, "\nFour most common organisms"),
         x = NULL, y = NULL) +
    theme_lab(10) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = if (n_years > 12) 90 else 0, vjust = 0.5, hjust = if (n_years > 12) 1 else 0.5))
}
ggsave(report_path("figures/isolates_per_year.png"), isolates_plot(FALSE), width = 8.5, height = 4.2, dpi = 150,
       bg = "white")
ggsave(file.path(gray_dir, paste0(run_name, "_isolates_per_year.png")), isolates_plot(TRUE), width = 8.5,
       height = 4.2, dpi = 110, bg = "white")

first_by_year <- first |>
  group_by(year) |>
  summarise(n = n(), n_cons = sum(mo == "B_STPHY_CONS"))

results <- list(params = params, n_isolates = nrow(isolates), n_first = nrow(first),
                first_by_year = first_by_year,
                abg_neg = abg_neg, abg_pos = abg_pos, by_period = by_period, trend = trend)
message("Wrote ", params$reports, "/tables, figures and antibiogram.pdf")
