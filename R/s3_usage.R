# Example Usage
input_data <- data.frame(
  Date = seq(as.Date("2020-01-01"), by = "week", length.out = 208),
  EQ_Volume = rnorm(208, mean = 100, sd = 20),
  Media_1 = rnorm(208, mean = 50, sd = 10),
  Media_2 = rnorm(208, mean = 30, sd = 5)
)

input_dates <- input_data$Date
input_dist_params <- list(
  Weibull = list(param1 = c(0.1, 5, 0.1), param2 = c(1, 70, 1)),
  Lognormal = list(param1 = c(1, 20, 1), param2 = c(10, 20, 1)),
  Logistic = list(param1 = c(1, 16, 1), param2 = c(1, 16, 1)),
  NegExp = list(param1 = c(0.1, 1, 0.1), param2 = c(NA, NA, NA)),
  ChiSquare = list(param1 = c(1, 70, 1), param2 = c(NA, NA, NA))
)

adstock_model <- AdstockModel(input_data, input_dates, list("linear", "diminishing", "log", "exponential"), input_dist_params)
adstock_model <- runAdstock(adstock_model)
write.csv(adstock_model$results, "adstock_new_iterations.csv")
