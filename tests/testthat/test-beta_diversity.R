library(PhyloIgSeq)

# ---- scree_plot() ----

test_that("scree_plot warns and returns NULL for NULL input", {
  expect_warning(result <- scree_plot(NULL), "No eigen values")
  expect_null(result)
})

test_that("scree_plot warns and returns NULL for empty input", {
  expect_warning(result <- scree_plot(numeric(0)), "No eigen values")
  expect_null(result)
})

test_that("scree_plot returns a ggplot for named eigenvalues", {
  eigen_values <- c(Axis.1 = 5, Axis.2 = 3, Axis.3 = 1, Axis.4 = 0.5)
  plt <- scree_plot(eigen_values)
  expect_s3_class(plt, "ggplot")
})

test_that("scree_plot computes percent variability against positive eigenvalues only", {
  eigen_values <- c(Axis.1 = 6, Axis.2 = 2, Axis.3 = 2)
  plt <- scree_plot(eigen_values)
  built <- ggplot2::layer_data(plt)
  # total_var = 10, so proportions are 60/20/20
  expect_equal(sort(built$y), c(20, 20, 60))
})

test_that("scree_plot excludes negative eigenvalues from the variance total but still plots them", {
  eigen_values <- c(Axis.1 = 8, Axis.2 = 2, Axis.3 = -1)
  plt <- scree_plot(eigen_values)
  built <- ggplot2::layer_data(plt)
  # total_var = 10 (only positive values count), Axis.3 still plotted at -10%
  expect_equal(built$y, c(80, 20, -10))
})

test_that("scree_plot truncates to max_nb_comp leading axes", {
  eigen_values <- c(Axis.1 = 5, Axis.2 = 4, Axis.3 = 3, Axis.4 = 2, Axis.5 = 1)
  plt <- scree_plot(eigen_values, max_nb_comp = 2)
  built <- ggplot2::layer_data(plt)
  expect_equal(nrow(built), 2)
})

test_that("scree_plot falls back to integer dim labels when eigen_values is unnamed", {
  eigen_values <- c(5, 3, 1)
  plt <- scree_plot(eigen_values)
  expect_equal(levels(plt$data$dim), c("1", "2", "3"))
})

test_that("scree_plot dim factor levels preserve input order", {
  eigen_values <- c(Axis.2 = 3, Axis.1 = 5, Axis.3 = 1)
  plt <- scree_plot(eigen_values)
  expect_equal(levels(plt$data$dim), c("Axis.2", "Axis.1", "Axis.3"))
})
