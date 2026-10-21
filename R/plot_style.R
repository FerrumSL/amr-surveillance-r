# Shared look for every figure, and the checks that keep it readable in print.
# Sourced by R/02_analyse.R and tests/testthat/test_style.R. Needs only farver (via ggplot2).

ink <- "#0b0b0b"
ink_2 <- "#52514e"
ink_3 <- "#a3a29d"  # "—" for too few isolates: deliberately faint
rule <- "#e4e3df"

# Colour from CIE L*, chroma and hue (polar Luv, the HCL model), so lightness is set exactly.
hcl_hex <- function(l, c, h) farver::encode_colour(cbind(h = h, c = c, l = l), from = "hcl")  # farver order: h, c, l

# % susceptible bins, light (most susceptible) to dark (least). One blue hue; lightness falls
# by 14 L* per step so every bin stays distinct in black-and-white print.
band_labels <- c(">= 90%", "80-89%", "60-79%", "40-59%", "20-39%", "< 20%")
band_breaks <- c(-Inf, 20, 40, 60, 80, 90, Inf)
band_fill <- setNames(hcl_hex(l = c(95, 81, 67, 53, 39, 25), c = c(12, 30, 45, 52, 50, 42), h = 250),
                      band_labels)
band_of <- function(pct) {
  factor(as.character(cut(pct, band_breaks, right = FALSE, labels = rev(band_labels))), band_labels)
}

# First isolates per year: four organisms and "Other", each a different hue AND lightness
# (L* 22-86 in steps of 16) so they also separate in grayscale. Other is the neutral middle.
category_fill <- c(hcl_hex(22, 45, 255), hcl_hex(38, 75, 25), hcl_hex(70, 55, 140),
                   hcl_hex(86, 70, 80))
other_fill <- hcl_hex(54, 0, 0)

lstar <- function(hex) farver::convert_colour(farver::decode_colour(hex), "rgb", "lab")[, "l"]

# WCAG 2 contrast ratio between two colours.
contrast <- function(a, b) {
  lum <- function(hex) {
    v <- farver::decode_colour(hex) / 255
    v <- ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
    drop(v %*% c(0.2126, 0.7152, 0.0722))
  }
  la <- lum(a)
  lb <- lum(b)
  (pmax(la, lb) + 0.05) / (pmin(la, lb) + 0.05)
}

# Black or white, whichever contrasts more with the fill.
text_on <- function(fill) ifelse(contrast(fill, "#000000") >= contrast(fill, "#ffffff"), "#000000", "#ffffff")

# Same lightness, no colour: what a black-and-white printer does.
to_gray <- function(hex) {
  hcl <- farver::decode_colour(hex, to = "hcl")
  farver::encode_colour(cbind(h = 0, c = 0, l = hcl[, "l"]), from = "hcl")
}

# Fails if neighbouring colours are closer than `min_gap` in L*, or if the text colour chosen
# for a fill has less than `min_contrast`. Returns the table it checked.
check_palette <- function(fill, min_gap = 12, min_contrast = 4.5, text = TRUE) {
  out <- data.frame(bin = names(fill), fill = unname(fill), L = round(lstar(fill), 1))
  out$gap_to_next <- c(round(abs(diff(out$L)), 1), NA)
  if (text) {
    out$text <- text_on(fill)
    out$contrast <- round(contrast(fill, out$text), 2)
  }
  bad_gap <- which(out$gap_to_next < min_gap)
  if (length(bad_gap) > 0) stop("L* gap below ", min_gap, " after: ", paste(out$bin[bad_gap], collapse = ", "))
  if (text && any(out$contrast < min_contrast)) {
    stop("Text contrast below ", min_contrast, ":1 on: ", paste(out$bin[out$contrast < min_contrast], collapse = ", "))
  }
  out
}

theme_lab <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = ink),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3)),
      plot.title.position = "plot", plot.caption.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = ink_2, margin = ggplot2::margin(b = 10), lineheight = 1.15),
      plot.caption = ggplot2::element_text(colour = ink_2, hjust = 0, size = ggplot2::rel(0.8),
                                           lineheight = 1.15, margin = ggplot2::margin(t = 10)),
      axis.text = ggplot2::element_text(colour = ink_2),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = rule, linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0, colour = ink),
      legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(14, 16, 12, 14))
}
