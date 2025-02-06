library(testthat)
library(dplyr)
library(data.table)

# Load the AdstockModel source file if needed
# source("adstock_model.R")

# Sample data for testing
sample_data <- data.frame(
  Date = seq(as.Date("2020-01-01"), by = "week", length.out = 10),
  EQ_Volume = rnorm(10, mean = 100, sd = 20),
  Media_1 = rnorm(10, mean = 50, sd = 10),
  Media_2 = rnorm(10, mean = 30, sd = 5)
)

sample_saturation_models <- list(
  linear = function(x) x,
  diminishing = function(x) x / (1 + x),
  log = function(x) log(1 + x),
  exponential = function(x) 1 - exp(-x)
)

sample_distribution_params <- list(
  Weibull = list(param1 = c(0.1, 5, 0.1), param2 = c(1, 70, 1)),
  Lognormal = list(param1 = c(1, 20, 1), param2 = c(10, 20, 1)),
  Logistic = list(param1 = c(1, 16, 1), param2 = c(1, 16, 1)),
  NegExp = list(param1 = c(0.1, 1, 0.1), param2 = c(NA, NA, NA)),
  ChiSquare = list(param1 = c(1, 70, 1), param2 = c(NA, NA, NA))
)

# Initialize Adstock Model
adstock_model <- AdstockModel(
  data = sample_data,
  dates = sample_data$Date,
  saturation_models = sample_saturation_models,
  distribution_params = sample_distribution_params
)

# Test initialization
test_that("Adstock model initializes correctly", {
  expect_true(is.list(adstock_model))
  expect_true("data" %in% names(adstock_model))
  expect_true("saturation_models" %in% names(adstock_model))
  expect_true("distribution_params" %in% names(adstock_model))
})

# Test running adstock
test_that("Adstock model runs without errors", {
  result <- runAdstock(adstock_model)
  expect_true(!is.null(result$results))
  expect_true(nrow(result$results) > 0)
})

# Test correlation computation
test_that("Correlation computation works correctly", {
  result <- runAdstock(adstock_model)
  correlations <- computeCorrelations(result)
  expect_true("Variable" %in% colnames(correlations))
  expect_true("Correl" %in% colnames(correlations))
  expect_true("Pval" %in% colnames(correlations))
})

# Test normalization function
test_that("Normalization function works correctly", {
  normalized_data <- normalize(c(1, 2, 3, 4, 5))
  expect_equal(sum(normalized_data), 1, tolerance = 1e-6)
})

# Test z-transformation function
test_that("Z-transformation standardizes correctly", {
  data <- c(1, 2, 3, 4, 5)
  transformed_data <- ztran(data)
  expect_equal(mean(transformed_data), 0, tolerance = 1e-6)
  expect_equal(sd(transformed_data), 1, tolerance = 1e-6)
})

# Run tests
test_file("test_adstock_model.R")
