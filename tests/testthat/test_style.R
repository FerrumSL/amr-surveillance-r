# Run from the project root: Rscript -e 'testthat::test_dir("tests/testthat")'
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "R/plot_style.R"), local = TRUE)

test_that("antibiogram bins stay apart in grayscale and their text is readable", {
  out <- check_palette(band_fill)
  expect_true(all(diff(out$L) < 0))  # lightness falls from >= 90% to < 20%
  expect_true(all(out$gap_to_next >= 12, na.rm = TRUE))
  expect_true(all(out$contrast >= 4.5))
  expect_true(all(abs(diff(lstar(to_gray(band_fill)))) >= 12))
  expect_error(check_palette(c(a = "#888888", b = "#8c8c8c")), "gap")
})

test_that("organism colours differ in lightness, not only hue", {
  all <- c(category_fill, other_fill)
  expect_true(all(diff(sort(lstar(all))) >= 12))
})

test_that("percentages fall into the documented bins", {
  expect_equal(as.character(band_of(c(100, 90, 89, 80, 79, 60, 59, 40, 39, 20, 19, 0))),
               rep(band_labels, each = 2))
})
